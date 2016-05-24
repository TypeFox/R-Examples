#******************************************************************************
#
# Copyright Tom Shott, ETA, PSU, 2013, 2014, 2015
# Use granted under BSD license terms
# R TFDEA Package
#
# NOTE: See utility.R for description of naming nomenclature, use of .l & .b suffix
# and other coding conventions
#
# Internal fast DEA function - no secondary, no lambdas, no duals, no slack maximization
# Does dea & sdea
#
#******************************************************************************
#
#library(lpSolveAPI)
# Load Utility Functions - included in Package
# source("utility.R")
# source("dea_common.R")

.dea_fast_internal <- function(x, y, rts="vrs", orientation="input",
                               super=TRUE,
                               round=FALSE, stdeff=FALSE,
                               index.K=NULL, index.T=NULL, debug=1){

  ###############################################################################################
  #
  # Setup results names, results output arrays, build technology set
  #
  # Note on xT, yT, and zT: these are the subet of DMU's that are the technology set compared with
  # Must use check index to convert logical to sparse index - TFDEA calls .sdea with logical index
  #
  ###############################################################################################
  # Inputs
  # index.K       <- .checkIndex(index.K, x, "index.K")     # DMU's to calc eff for
  # index.T       <- .checkIndex(index.T, x, "index.T")     # Technology Set

  nd            <- nrow(x)                  # number of units, firms, DMUs
  names.dmu     <- rownames(x)
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

  # Outputs
  dmu.status1    <- array(NA, c(nd),        list(dmu=names.dmu))
  efficiency    <- array(NA, c(nd),         list(dmu=names.dmu))

  ###############################################################################################
  #
  # Setup constant parts of linear equation that don't change per K. Each DMU uses same values
  #   other parts of equation need to change before we solve for each K DMU
  #
  ###############################################################################################
  # Step 1 - Create model
  lpm <- .dea_create_model(orientation=orientation, rts=rts, nT=nT, nv=nv,
                           vars.super=TRUE, vars.cook=FALSE, vars.slack=FALSE, vars.phase2=FALSE,
                           rname=c(names.in, names.out), cname=names.dmu[index.T], debug=debug)

  # Step 2 - Setup control values
  lpm           <- .dea_set_ctl(lpm, debug=debug)
  epsilon       <- lpm$epsilon

  # Step 3- Calculate and set objfun
  obj.phase1    <- .dea_obj_phase1(lpm, vars.cook=FALSE, debug=debug)
  set.objfn(lpm$model, obj.phase1)

  # Steps 4 & 5 - setup RTS & Technology Set of model
  lpm           <- .dea_set_rts(lpm, debug=debug)
  lpm           <- .dea_set_technology(lpm, xT=xT, yT=yT, vars.slack=FALSE, debug=debug)

  # Step 5 - Precalculate values for inner loop for eff col, rhs col and super row
  values.eff    <- .dea_values_eff(lpm, x=x, y=y, debug=debug)
  values.rhs    <- .dea_values_rhs(lpm, x=x, y=y, debug=debug)
  values.super  <- .dea_values_super(lpm, nd=nd, index.K=index.K, index.T=index.T,
                                     super=super, debug=debug)

  set.basis(lpm$model, default=TRUE)

  if (debug >= 2) lpDebugFile(lpm$model, "Done constant part eqn")

  ###############################################################################################
  #
  # For Each DMU K change DMU specific values and then solve equation & save results
  # Loop for each DMU - Solve one linear equation and extract eff for each DMU
  #
  ###############################################################################################

  for (k in index.K)  {
    # setup eff col, rhs col & super row for each K
    # setting col clears obj fun - must inc. 0 value for objfun in lpm.col.eff to not clear
    set.column(lpm$model, lpm$col.eff,      values.eff[k,],       c(0:nv))   # keep onjfun
    set.rhs(lpm$model,                      values.rhs[k,],       c(1:nv))
    set.row(lpm$model,    lpm$row.super,    values.super[k,],     c(1:lpm$col.n))

    #     if (debug >= 4)
    #       lpDebugFile(lpm$model, paste0("Final constraints added k= ",k)) # print for every K

    # Solve LP model - Phase 1
    dmu.status1[k] <- solve.lpExtPtr(lpm$model)
    efficiency[k] <- get.variables(lpm$model)[lpm$col.eff]
  }

  ###############################################################################################
  #
  # Cleanup & return Results
  #
  ###############################################################################################
  # NOTE: shape & names of ifelse test must be same as var to retain names

  for(i in which(dmu.status1 != 0)){
    .print_status_msg(dmu.status1[i], 1, i, debug=debug)
  }


  # Cleanup Eff values based on solver status
  efficiency <- ifelse(dmu.status1==0,
                       efficiency,
                       ifelse( dmu.status1 == 2 | dmu.status1 == 3,
                        ifelse(array(lpm$orientation.in, c(nd), list(dmu=names.dmu)), Inf, -Inf),
                        NA))

  # WARN: Rounding turned off by default. To match other DEA packages, turn rounding on,
  # but does hide some numerical results.
  # Round efficiencies that are within epilson to 1 or zero
  efficiency[abs(efficiency-1) < epsilon & round] <- 1
  efficiency[abs(efficiency)   < epsilon & round] <- 0

  efficiency  <- ifelse( array(stdeff && !lpm$orientation.in, c(nd), list(dmu=names.dmu)),
                       1 / efficiency, efficiency)

  results.l   <- list(eff=efficiency, dmu.status1=dmu.status1)

  return(results.l)
}

#<new Page>
###############################################################################################
#
# Wrapper with input checking for internal fast DEA / SDEA function
#
###############################################################################################
#
.dea_fast <- function(x, y, rts="vrs", orientation="input", super=FALSE,
                      round=FALSE, stdeff=FALSE,
                      index.K=NULL, index.T=NULL, debug=1) {

  rts         <- .checkOption(rts,           "rts",         options.rts.l)
  orientation <- .checkOption(orientation,   "orientation", options.orientation.l)

  super       <- .checkOption(super,         "super",       TRUE)

  round       <- .checkOption(round,         "round",       TRUE)
  stdeff      <- .checkOption(stdeff,        "stdeff",      TRUE)
  debug       <- .checkOption(debug,         "debug",       0)

  # Check that x & y are legal inputs & convert to standard dimension arrays
  x <- .checkData(x, "x")
  y <- .checkData(y, "y")

  # .checkDataGood now includes this check
  #   if (nrow(x) != nrow(y))
  #     stop("Number of DMU's in inputs != number of DMU's in outputs", call. = FALSE)
  .checkDataGood(x, y, debug = max(0, debug-1))

  # There functions also convert logical index to sparse index
  index.K       <- .checkIndex(index.K, x, "index.K")     # DMU's to calc eff for
  index.T       <- .checkIndex(index.T, x, "index.T")     # Technology Set

  results <- .dea_fast_internal(x, y, rts=rts, orientation=orientation,
                       super=super, round=round, stdeff=stdeff,
                       index.K=index.K, index.T=index.T,
                       debug=debug)

  return(results)
}
