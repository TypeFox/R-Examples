#******************************************************************************
#
# Copyright Tom Shott, ETA, PSU, 2013, 2014
# Use granted under BSD license terms
#
# R TFDEA Package
#
# NOTE: See utility.R for description of naming nomenclature, use of .l & .b suffix
# and other coding conventions
#
# Internal DEA function - not public
#
#******************************************************************************
#
# Load Utility Functions
# source("utility.R")
# source("dea_common.R")


.dea_internal <- function(x, y, rts="vrs", orientation="input",
                          slack=TRUE, dual=FALSE,
                          second.b=FALSE, zn=NULL,
                          round=FALSE, stdeff=FALSE,
                          index.K=NULL, index.T=NULL, debug=1){

  ###############################################################################################
  #
  # Setup results names, results output arrays, build technology set
  #
  # Note on xT, yT, and zT: these are the subet of DMU's that are the technology set compared with
  # Must call with sparse index - .dea converts sparse to logical
  # No input type checking, all done in .dea
  #
  ###############################################################################################
  # Inputs

  nd            <- nrow(x)                  # number of units, firms, DMUs
  names.dmu     <- rownames(x)

  index.K       <- ifelse(rep.int(is.null(index.K), nd), c(1:nd), index.K)
  index.T       <- ifelse(rep.int(is.null(index.T), nd), c(1:nd), index.T)
  nT            <- length(index.T)          # Can't use nrow, DMU length 1 has no rows

  nx            <- ncol(x)                  # number of inputs
  names.in      <- colnames(x)

  # WARNING: Must use drop option to prevent dimensions from being reduced
  xT            <- x[index.T, ,drop=FALSE]  # Technology Set - inputs

  ny            <- ncol(y)                  # number of outputs
  names.out     <- colnames(y)
  yT            <- y[index.T, ,drop=FALSE]  # Technology Set - Outputs

  nv            <- nx + ny

  zn            <- ifelse(rep.int(is.null(zn), nd),    0,    zn)
  zT            <- zn[index.T, drop=FALSE]   # Technology Set - secondary objective

  efficiency    <- array(NA, c(nd),       list(dmu=names.dmu))
  dmu.status1   <- array(NA, c(nd),       list(dmu=names.dmu))
  dmu.status2   <- array(NA, c(nd),       list(dmu=names.dmu))

  lambda        <- array(NA, c(nd, nd),   list(dmu=names.dmu, dmu2=names.dmu))

  if (dual)
    weight.xy   <- array(NA, c(nd, nv+1), list(dmu=names.dmu, "xy"= c("w", names.in, names.out)))

  if(slack)
    slack.xy    <- array(NA, c(nd, nv),   list(dmu=names.dmu, "xy"= c(names.in, names.out)))


  ###############################################################################################
  #
  # Setup constant parts of linear equation that don't change per K
  #
  ###############################################################################################
  # Step 1 - Create model

  phase2.run.b  <- slack || second.b
  vars.slack.b  <- slack

  lpm <- .dea_create_model(orientation=orientation, rts=rts, nT=nT, nv=nv,
                           vars.slack=vars.slack.b, vars.phase2=phase2.run.b,
                           rname=c(names.in, names.out), cname=names.dmu[index.T],
                           debug=debug)

  # Step 2 - Setup control values
  lpm           <- .dea_set_ctl(lpm, debug=debug)
  epsilon       <- lpm$epsilon

  # Step 3- Calculate and set objfuns
  obj.phase1    <- .dea_obj_phase1(lpm, debug=debug)
  set.objfn(lpm$model, obj.phase1)

  obj.phase2    <- .dea_obj_phase2(lpm, zT=zT, vars.slack=vars.slack.b, debug=debug)

  # Steps 4 & 5 - setup RTS & Technology Set of model
  lpm           <- .dea_set_rts(lpm, debug=debug)
  lpm           <- .dea_set_technology(lpm, xT=xT, yT=yT, vars.slack=vars.slack.b, debug=debug)

  # Step 5 - Precalculate values for inner loop for eff col, rhs col and super row
  values.eff    <- .dea_values_eff(lpm, x=x, y=y, debug=debug)
  values.rhs    <- .dea_values_rhs(lpm, x=x, y=y, debug=debug)

  ###############################################################################################
  #
  # For Each K setup values and solve equation
  # Loop for each DMU - Solve one linear equation and extract eff for each DMU
  #
  ###############################################################################################
  #
  for (k in index.K){
    #
    # Setup Variable part LP model that changes for each K
    #

    # setting col clears obj fun - must reset objfun
    set.column(lpm$model, lpm$col.eff,      values.eff[k,],         c(0:nv))  # Set objfun also
    set.rhs(lpm$model,                      values.rhs[k,],         c(1:nv))

    tmp.status          <- solve.lpExtPtr(lpm$model)
    dmu.status1[k]      <- tmp.status

    tmp.var             <- get.variables(lpm$model)[1:lpm$col.n]
    phase1.eff          <- tmp.var[lpm$col.eff]
    efficiency[k]       <- phase1.eff

    if (dual){
      weight.xy[k,]     <- ifelse( rep.int(lpm$orientation.in, nv+1),
                                   get.dual.solution(lpm$model), -get.dual.solution(lpm$model))
    }

    #
    # Phase 2 - slack maximization or secondary objective
    #
    if (tmp.status == 0 && phase2.run.b){

      phase1.eff <- ifelse( abs(phase1.eff-1) < epsilon, 1, phase1.eff)     # Round to 1
      phase1.eff <- ifelse( abs(phase1.eff)   < epsilon, 0, phase1.eff)     # Round to 0

      # Add constraint that efficiency must equal efficiency found in Phase 1
      set.row(lpm$model,                (lpm$row.eff), 1, c(lpm$col.eff))   # Lock eff
      set.constr.type(lpm$model, "=",  c(lpm$row.eff) )                     # Make constraint =
      set.rhs(lpm$model, phase1.eff,   c(lpm$row.eff) )                     # Make eff = phase 1

      set.objfn(lpm$model, obj.phase2)

      if (k == 1 && debug >= 2) lpDebugFile(lpm$model, "Phase 2 - Pre slack maximizations")

      tmp.status      <- solve.lpExtPtr(lpm$model)
      dmu.status2[k]  <- tmp.status
      tmp.var         <- get.variables(lpm$model)[1:lpm$col.n]

      if (slack){
        col.slack     <- lpm$col.slack[1]
        slack.xy[k,]  <- tmp.var[ col.slack   : (col.slack + nv - 1)]
      }

      # Cleanup from phase 2
      set.constr.type(lpm$model, 0, c(lpm$row.eff))     # Make slack constraint free
      set.objfn(lpm$model, obj.phase1)                  # Need to clear out slack  values

    } # End Phase 2

    # Save last lambda value calculated
    lambda[k, index.T]  <- tmp.var[1:nT]

  } # End of DMU K Loop


  ###############################################################################################
  #
  # Cleanup & return Results
  #
  ###############################################################################################

  for(k in which(dmu.status1 != 0)){
    .print_status_msg(dmu.status1[k], 1, k, debug=debug)
  }

  for(i in which(!is.na(dmu.status2) & dmu.status2 != 0)){
    .print_status_msg(dmu.status2[i], 2, i, debug=debug)
  }

  dmu.error.b <- (dmu.status1 != 0) | (!is.na(dmu.status2) & dmu.status2 != 0)

  efficiency <- ifelse(dmu.status1 == 0,
                       efficiency,
                       ifelse( dmu.status1 == 2 | dmu.status1 == 3,
                               ifelse(lpm$orientation.in, Inf, -Inf),
                               NA))


  # WARN: Rounding turned off by default. To match other DEA packages, turn rounding on,
  # but does hide some numerical results.

  # Round efficiencies that are within epilson to 1 or zero to 1 or zero.
  efficiency[abs(efficiency-1) < epsilon & rep.int(round, nd)] <- 1
  efficiency[abs(efficiency)   < epsilon & rep.int(round, nd)] <- 0

  efficiency  <- ifelse( array(stdeff && !lpm$orientation.in, c(nd), list(dmu=names.dmu)),
                         1 / efficiency, efficiency)

  lambda[dmu.error.b, ] <- NA
  results.l   <- list(eff=efficiency, status1=dmu.status1, status2=dmu.status2, lambda=lambda)

  if (dual){
    weight.xy[dmu.status1 != 0] <- NA

    weight.w    <- ifelse(array(lpm$orientation.in, c(nd)), weight.xy[ 1 ], -weight.xy[ 1 ])
    weight.x    <- weight.xy[ , c(     2  : (nx+1)), drop=FALSE ]
    weight.y    <- weight.xy[ , c( (nx+2) : (nv+1)), drop=FALSE ]

    dimnames(weight.w)  <- list(dmu=names.dmu)
    dimnames(weight.x)  <- list(dmu=names.dmu, vx=names.in)
    dimnames(weight.y)  <- list(dmu=names.dmu, uy=names.out)

    results.l           <- c(results.l, list(vx=weight.x, uy=weight.y, w=weight.w))
  }


  if(slack){
    slack.xy[ !is.na(dmu.status2) & dmu.status2 != 0 ] <- NA

    slack.x <- slack.xy[ , c( 1      : nx), drop=FALSE]
    slack.y <- slack.xy[ , c( (nx+1) : nv), drop=FALSE]

    dimnames(slack.x)   <- list(dmu=names.dmu, sx=names.in)
    dimnames(slack.y)   <- list(dmu=names.dmu, sy=names.out)

    results.l           <- c(results.l, list(sx=slack.x, sy=slack.y))
  }

  return(results.l)
}


# End dea function


#<New Page>
#
# Wrapper with input checking for internal DEA functions and use has options external DEA function
# does not include
#
.dea <- function(x, y, rts="vrs", orientation="input",
                          slack=TRUE, dual=FALSE,
                          second="none", z=0,
                          round=FALSE, stdeff=FALSE,
                          index.K=NULL, index.T=NULL,
                          debug=1) {

  rts         <- .checkOption(rts,           "rts",         options.rts.l)
  orientation <- .checkOption(orientation,   "orientation", options.orientation.l)

  slack       <- .checkOption(slack,         "slack",       TRUE)
  dual        <- .checkOption(dual,          "dual",        TRUE)

  second      <- .checkOption(second,        "second",      options.second.l)

  round       <- .checkOption(round,         "round",       TRUE)
  stdeff      <- .checkOption(stdeff,        "stdeff",      TRUE)
  debug       <- .checkOption(debug,         "debug",       0)

  # Check for legal input values
  #
  # Check that x & y are legal inputs & convert to standard dimension arrays
  x <- .checkData(x, "x")
  y <- .checkData(y, "y")

  # .checkDataGood now includes this check
  #   if (nrow(x) != nrow(y))
  #     stop("Number of DMU's in inputs != number of DMU's in outputs", call. = FALSE)
  .checkDataGood(x, y, debug = max(0 , debug-1))


  # Must use check index to convery logical to sparse index
  # TFDEA calls .dea with logical index
  index.K       <- .checkIndex(index.K, x, "index.K")     # DMU's to calc eff for
  index.T       <- .checkIndex(index.T, x, "index.T")     # Technology Set

  zn <- rep.int(0, nrow(x))
  second.b <- (second == "min" || second == "max")

  # Check secondary optimization parms
  if (second.b){
    z <- .checkVector(z, "z")
    if (nrow(x) != length(z))
      stop("secondary data size (rows) must match number of DMU's", call. = FALSE)

    if (slack)
      stop("Can not use both second and slack option; add option slack=FALSE", call. = FALSE)

    zn <- ifelse( rep(second == "min", nrow(x)),  -z,  z)
  }

  results <- .dea_internal(x, y, rts=rts, orientation=orientation,
                           slack=slack, dual=dual,
                           second.b=second.b, zn=zn,
                           round=round, stdeff=stdeff,
                           index.K=index.K, index.T=index.T,
                           debug=debug)

  return(results)

}

