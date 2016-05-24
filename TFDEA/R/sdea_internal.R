#******************************************************************************
#
# Copyright Tom Shott, ETA, PSU, 2013, 2014, 2015
# Use granted under BSD license terms
#
# R TFDEA Package
#
# NOTE: See utility.R for description of naming nomenclature, use of .l & .b suffix
# and other coding conventions
#
# Internal SDEA function
#
# Update version that does super efficiency calculations using Wade Cook Alog.
# ToDo: Add citation Wade Cook paper
#
#******************************************************************************
#
# library(lpSolveAPI)
# Load Utility Functions
# source("utility.R")
# source("dea_common.R")

#
# Internal function - not public
#
.sdea_internal <- function(x, y, rts="vrs", orientation="input", slack=TRUE, dual=FALSE,
                  cook=FALSE, second.b=FALSE, zn=NULL,
                  round=FALSE, stdeff=FALSE,
                  index.K=NULL, index.T=NULL, debug=1){

  ###############################################################################################
  #
  # Setup results names, results output arrays, build technology set
  #
  # Note on xT, yT, and zT: these are the subset of DMU's that are the technology set compared with
  # Must use check index to convert logical to sparse index - TFDEA calls .sdea with logical index
  #
  ###############################################################################################
  # Inputs

  nd            <- nrow(x)                  # number of units, firms, DMUs
  names.dmu     <- rownames(x)
  index.K       <- ifelse(rep.int(is.null(index.K), nd), c(1:nd), index.K)
  index.T       <- ifelse(rep.int(is.null(index.T), nd), c(1:nd), index.T)


  nT            <- length(index.T)          # Can't use nrow, DMU length 1 has no rows
  nK            <- length(index.K)

  nx            <- ncol(x)                  # number of inputs
  names.in      <- colnames(x)

  # WARNING: Must use drop option to prevent dimensions from being reduced
  xT            <- x[index.T, ,drop=FALSE]  # Technology Set - inputs

  ny            <- ncol(y)                  # number of outputs
  names.out     <- colnames(y)
  yT            <- y[index.T, ,drop=FALSE]  # Technology Set - Outputs

  nv            <- nx + ny                  # Number of vars

  zn            <- ifelse(rep.int(is.null(zn), nd),    0,    zn)
  zT            <- zn[index.T, drop=FALSE]   # Technology Set - secondary objective


  # Outputs
  dmu.status1   <- array(NA, c(nd),     list(dmu=names.dmu))
  dmu.status2   <- array(NA, c(nd),     list(dmu=names.dmu))
  efficiency    <- array(NA, c(nd),     list(dmu=names.dmu))  # efficiency

  se.excess     <- array(NA, c(nd),     list(dmu=names.dmu))  # Cook Super efficiency Vars
  se.eff        <- array(NA, c(nd),     list(dmu=names.dmu))

  lambda        <- array(NA, c(nd, nd), list(dmu=names.dmu, dmu2=names.dmu))


  if (dual){
    weight.xy   <- array(NA, c(nd, (nv+1)),
                         list(dmu=names.dmu, "xy"= c("w", names.in, names.out)))
  }

  if (slack){
    slack.xy    <- array(NA, c(nd, nv),
                         list(dmu=names.dmu, "xy"= c(names.in, names.out)))
  }

  ###############################################################################################
  #
  # Setup constant parts of linear equation that don't change per K
  #
  ###############################################################################################

  phase2.run.b  <- slack || second.b
  vars.slack.b  <- slack

  # Step 1 - Create model
  lpm <- .dea_create_model(orientation=orientation, rts=rts, nT=nT, nv=nv,
                           vars.super=TRUE, vars.cook=TRUE, vars.slack=vars.slack.b,
                           vars.phase2=phase2.run.b,
                           rname=c(names.in, names.out), cname=names.dmu[index.T], debug=debug)

  # Step 2 - Setup control values
  lpm         <- .dea_set_ctl(lpm, debug=debug)
  epsilon     <- lpm$epsilon

  # Step 3- Calculate and set objfn
  obj.phase1  <- .dea_obj_phase1(lpm, vars.cook=TRUE, debug=debug)
  set.objfn(lpm$model, obj.phase1)

  obj.phase2  <- .dea_obj_phase2(lpm, zT=zT, vars.slack=vars.slack.b, debug=debug)

  # Steps 4 & 5 - setup RTS & Technology Set of model
  lpm         <- .dea_set_rts(lpm, debug=debug)
  lpm         <- .dea_set_technology(lpm, xT=xT, yT=yT, vars.slack=vars.slack.b, debug=debug)

  # If no phase2 - still need to lock phase1 for non Cook
  set.row(lpm$model, lpm$row.exc, as.numeric(!cook), c(lpm$col.exc))
  #set.constr.type(lpm$model, c(0),      c(lpm$row.eff) )             # Make eff free constraint

  # Step 5 - Precalculate values for inner loop for eff col, rhs col and super row
  values.eff  <- .dea_values_eff(lpm, x=x, y=y, debug=debug)
  values.exc  <- .dea_values_exc(lpm, x=x, y=y, debug=debug)
  values.rhs  <- .dea_values_rhs(lpm, x=x, y=y, cook=cook, debug=debug)
  values.super <- .dea_values_super(lpm, nd=nd, index.K=index.K, index.T=index.T,
                                     super=TRUE, debug=debug)

  ###############################################################################################
  #
  # For Each K setup values and solve equation
  # Loop for each DMU - Solve one linear equation and extract eff for each DMU
  #
  ###############################################################################################

  for (k in index.K)  {
    #
    # Setup Variable part LP model that changes for each K
    #

    # setting col clears obj fun - must reset objfn
    set.column(lpm$model, lpm$col.eff,      values.eff[k,],         c(0:nv))  # Set objfn also
    set.column(lpm$model, lpm$col.exc,      values.exc[k,],         c(0:nv))
    set.rhs(lpm$model,                      values.rhs[k,],         c(1:nv))
    set.row(lpm$model,    lpm$row.super,    values.super[k,],       c(1:lpm$col.n))

    # Cleared by col set - must be re-locked
    set.row(lpm$model, lpm$row.exc, as.numeric(!cook), c(lpm$col.exc)) # Lock cook for normal

    if ((debug >= 3 && k==1 ) || debug >= 4)
      lpDebugFile(lpm$model, paste0("Post-set obj constraints added k= ",k)) # print for every K

    #
    # Solve LP model - Phase 1
    #
    tmp.status      <- solve.lpExtPtr(lpm$model)
    dmu.status1[k]  <- tmp.status

    tmp.var         <- get.variables(lpm$model)[1:lpm$col.n]
    tmp.eff         <- tmp.var[lpm$col.eff]
    tmp.exc         <- tmp.var[lpm$col.exc]

    # NOTE: Cook paper uses a min eqn for input orientation and min eqn for output orientation.
    # This code implements input orientation as min and output orientation as max - thus
    # we need to invert the sign of excess for output orientation

    phase1.eff <- ifelse(!cook,
                         tmp.eff,                                                 # Std SDEA
                         ifelse(! (is.finite(tmp.exc) && abs(tmp.exc) > epsilon),
                                1 + tmp.eff,                                      # Cook No excess
                                ifelse(lpm$orientation.in,
                                       1 + tmp.eff + 1 / ( 1 - tmp.exc),          # Cook Yes excess
                                       1 / (1 + tmp.exc + 1 / (1 + tmp.eff)))))

    se.eff[k]           <- ifelse(lpm$orientation.in, tmp.eff, -tmp.eff)
    se.excess[k]        <- tmp.exc
    efficiency[k]       <- phase1.eff

    # If dual option, save dual values
    if (dual){
      weight.xy[k,]     <- ifelse( rep.int(lpm$orientation.in, nv+1),
                                   get.dual.solution(lpm$model), -get.dual.solution(lpm$model))
    }

    #
    # Phase 2 - slack maximization or secondary objective
    #

    # Need to lock results for (efficiency, excess) found in phase 1
    # Need to use actual values from model - not phase1 reported efficiency
    if (phase2.run.b && tmp.status == 0){

      tmp.eff <- ifelse( abs(tmp.eff-1) < epsilon, 1, tmp.eff)       # Round to 1
      tmp.eff <- ifelse( abs(tmp.eff)   < epsilon, 0, tmp.eff)       # Round to 0

      # Add constraint that efficiency & excess must equal Phase 1 values. Some cleared by
      # the col set in setting K in Phase 1
      set.constr.type(lpm$model, "=",         c(lpm$row.eff))    # = constraint
      set.row(lpm$model,                        lpm$row.eff, 1, lpm$col.eff)  # Lock eff value

      set.row(lpm$model,                        lpm$row.exc, 1, lpm$col.exc)  # Put 1 in for exc

      set.rhs(lpm$model,  c(tmp.eff,tmp.exc), c(lpm$row.eff, lpm$row.exc))    # Make  = phase 1

      set.objfn(lpm$model, obj.phase2)
      if (k == 1 && debug >= 2) lpDebugFile(lpm$model, "Phase 2 - Pre slack maximizations")

      #
      # Solve Model - Phase 2
      #
      tmp.status      <- solve.lpExtPtr(lpm$model)
      dmu.status2[k]  <- tmp.status
      tmp.var         <- get.variables(lpm$model)[1:lpm$col.n]

      if (slack){
        col.slack <- lpm$col.slack[1]
        slack.xy[k,] <- tmp.var[ col.slack         : (col.slack + nv - 1)]
      }

      # Undo phase 2 lock vars
      set.constr.type(lpm$model, c(0),      c(lpm$row.eff) )    # Make free constraint
      set.rhs(lpm$model,  0,                c(lpm$row.exc) )    # Put back to zero
      set.objfn(lpm$model, obj.phase1)                          # Reset objfn

    } # End Phase 2

    # Grab last version of lambda
    lambda[k, index.T]  <- tmp.var[1:nT]

  } # End of DMU Loop

  ###############################################################################################
  #
  # Cleanup & return Results
  #
  ###############################################################################################

  for(k in which(dmu.status1 != 0)){
    .print_status_msg(dmu.status1[k], 1, k, debug=debug)
  }
  for(k in which(dmu.status2 != 0 || dmu.status2 != NA)){
    .print_status_msg(dmu.status2[k], 2, k, debug=debug)
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

  efficiency    <- ifelse( array(stdeff && !lpm$orientation.in, c(nd), list(dmu=names.dmu)),
                       1 / efficiency, efficiency)

  results.l     <- list(eff=efficiency, status1=dmu.status1, status2=dmu.status2)

  if(cook){
    se.excess[dmu.status1 != 0] <- 0
    results.l   <- c(results.l, list(se.eff=se.eff, se.excess=se.excess))
  }

  lambda[dmu.error.b, ] <- NA
  diag(lambda)  <- 0                                 # Match Benchmarking - sdea self lambda zero
  results.l     <- c(results.l, list(lambda=lambda))

  if (dual){
    weight.xy[(dmu.status1 != 0)] <- NA

    weight.w    <- ifelse(array(lpm$orientation.in, c(nd)),  weight.xy[ 1 ],  -weight.xy[ 1 ])
    weight.x    <- weight.xy[ , c(     2  : (nx+1)), drop=FALSE ]
    weight.y    <- weight.xy[ , c( (nx+2) : (nv+1)), drop=FALSE ]

    dimnames(weight.w) <- list(dmu=names.dmu)
    dimnames(weight.x) <- list(dmu=names.dmu, vx=names.in)
    dimnames(weight.y) <- list(dmu=names.dmu, uy=names.out)

    results.l       <- c(results.l, list(vx=weight.x, uy=weight.y, w=weight.w))
  }

  if (slack){
    slack.xy[!is.na(dmu.status2) & dmu.status2 != 0] <- NA

    slack.x <- slack.xy[ , c( 1      : nx), drop=FALSE]
    slack.y <- slack.xy[ , c( (nx+1) : nv), drop=FALSE]

    dimnames(slack.x) <- list(dmu=names.dmu, sx=names.in)
    dimnames(slack.y) <- list(dmu=names.dmu, sy=names.out)

    results.l   <- c(results.l, list(sx=slack.x, sy=slack.y))
  }

  return(results.l)
}


#<New Page>
###############################################################################################
#
# Wrapper with input checking for internal SDEA function
#
###############################################################################################
#
.sdea <- function(x, y, rts="vrs", orientation="input",
                       slack=TRUE, dual=FALSE,
                       second="none", z=0,
                       cook=FALSE,
                       round=FALSE, stdeff=FALSE, fast=FALSE,
                       index.K=NULL, index.T=NULL, debug=1) {

  rts         <- .checkOption(rts,           "rts",         options.rts.l)
  orientation <- .checkOption(orientation,   "orientation", options.orientation.l)

  slack       <- .checkOption(slack,         "slack",       TRUE)
  dual        <- .checkOption(dual,          "dual",        TRUE)

  cook        <- .checkOption(cook,          "cook",        TRUE)

  second      <- .checkOption(second,        "second",      options.second.l)

  round       <- .checkOption(round,         "round",       TRUE)
  fast        <- .checkOption(fast,          "fast",        TRUE)
  stdeff      <- .checkOption(stdeff,        "stdeff",      TRUE)
  debug       <- .checkOption(debug,         "debug",       0)

  # Check that x & y are legal inputs & convert to standard dimension arrays
  x <- .checkData(x, "x")
  y <- .checkData(y, "y")

  # .checkDataGood now includes this check
  #   if (nrow(x) != nrow(y))
  #     stop("Number of DMU's in inputs != number of DMU's in outputs", call. = FALSE)
  .checkDataGood(x, y, debug = max(0 , debug-1))

  index.K       <- .checkIndex(index.K, x, "index.K")     # DMU's to calc eff for
  index.T       <- .checkIndex(index.T, x, "index.T")     # Technology Set

  zn <- rep.int(0, nrow(x))
  second.b <- (second == "min" || second == "max")

  # Check secondary optimization parms
  if (second.b){
    z <- .checkVector(z,"z")
    if (nrow(x) != length(z))
      stop("secondary data size (rows) must match number of DMU's", call. = FALSE)

    zn <- ifelse( rep.int(second == "min", nrow(x)),  -z,  z)
    if (slack)
      stop("Can not use both second and slack option; set slack=FALSE", call. = FALSE)
  }

  results <- .sdea_internal(x, y, rts=rts, orientation=orientation, slack=slack, dual=dual,
                            cook=cook, second.b=second.b, zn=zn,
                            round=round, stdeff=stdeff,
                            index.K=index.K, index.T=index.T, debug=debug)

  return(results)
}
