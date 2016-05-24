#
# Common DEA functions
#
###########################################################################################
#
# Create a dea model object - Step 1
# Size of model depends on if using slack vars, phase2 and cook option
#
###########################################################################################

.dea_create_model <- function(orientation, rts, nT, nv,
                              vars.super=FALSE, vars.cook=FALSE, vars.slack=FALSE, vars.phase2=FALSE,
                              rname, cname, debug=1){

  lpm <- list()

  lpm$orientation.in  <- orientation == "input"
  lpm$rts             <- rts

  lpm$nT              <- nT

  lpm$vars.super      <- vars.super
  lpm$vars.cook       <- vars.cook
  lpm$vars.slack      <- vars.slack
  lpm$vars.phase2     <- vars.phase2


  # Setup cols - - don't change order - some code uses adjacent cols
  lpm$col.n           <- nT
  lpm$col.T           <- c(1:nT)

  if (vars.cook){
    lpm$col.n         <- lpm$col.n + 1
    lpm$col.exc       <- lpm$col.n
    cname             <- c(cname, "exc")
  } else {
    lpm$col.exc       <- 1
  }

  lpm$col.n           <- lpm$col.n + 1
  lpm$col.eff         <- lpm$col.n
  cname               <- c(cname, "eff")

  if (vars.slack){
    lpm$col.slack     <- c( (lpm$col.n+1) : (lpm$col.n + nv) )
    lpm$col.n         <- lpm$col.n + nv
    cname             <- c(cname, paste("S", rname, sep="_"))
  }

  # Setup rows - don't change order - some code uses adjacent rows
  lpm$row.n           <- nv + 1
  lpm$row.rts         <- lpm$row.n
  rname               <- c(rname, "rts")
  lpm$row.lock        <- NULL

  if (vars.super){
    lpm$row.n         <- lpm$row.n + 1
    lpm$row.super     <- lpm$row.n
    rname             <- c(rname, "super")
    lpm$row.lock      <- c(lpm$row.n, lpm$row.super)
  }

  if(vars.phase2){
    lpm$row.n         <- lpm$row.n + 1
    lpm$row.eff       <- lpm$row.n
    lpm$row.lock      <- c(lpm$row.lock, lpm$row.eff)
    rname             <- c(rname, "eff")
  }

  if (vars.cook){                                         # If using cook, need to lock in Phase1
    lpm$row.n       <- lpm$row.n + 1
    lpm$row.exc     <- lpm$row.n
    lpm$row.lock    <- c(lpm$row.lock, lpm$row.exc)
    rname           <- c(rname, "exc")
  }


  lpm$model <- make.lp(lpm$row.n, lpm$col.n, verbose = "neutral")

  lpm.name <- paste0(getSrcFilename(DEA), orientation, rts, tfdea.id, sep="-")
  name.lp(lpm$model)

  dimnames(lpm$model)   <- list(rname, cname)

  # Set zero lower bound, except for Eff which is -Inf
  bounds                <- rep.int(0.0, lpm$col.n)
  bounds[lpm$col.eff]   <- -Inf
  set.bounds(lpm$model, lower = bounds)

  # setup standard constraints for non-Var constraints (eff, exc, super)
  if (length(lpm$row.lock) >= 1){
    set.constr.type(lpm$model,      rep.int("=", length(lpm$row.lock) ),  lpm$row.lock)
    set.rhs(lpm$model,              rep.int(0,   length(lpm$row.lock) ),  lpm$row.lock)
  }

  if (debug >= 2) lpDebugFile(lpm$model, "Model created")

  return(lpm)
}

###########################################################################################
#
# lpm setup control options - Step 2
#
###########################################################################################

.dea_set_ctl <-function (lpm, debug=1){

  lpm.sense       <- ifelse(lpm$orientation.in,     "min",      "max")
  lpm.ctl.opts    <- list(sense=lpm.sense, timeout=60)

  verbose.l       <- c("neutral", "critical", "severe", "important", "normal", "detailed")

   if (debug >= 1){
     debug <- min(debug, 6)
     lpm.ctl.opts <- c(lpm.ctl.opts, list(verbose = verbose.l[debug]))
   }

  lpm.ctl.opts <- c(lpm.ctl.opts, list(scaling = c("extreme", "quadratic")))
  lpm.ctl.opts <- c(lpm.ctl.opts, list(pivoting = c("devex", "loopalternate")))

  if (debug >= 3){
    cat("\nControl Options:")
    print(lpm.ctl.opts)
  }
  do.call( lp.control, c( list( lprec = lpm$model ), lpm.ctl.opts ) )
  lpm$ctl <- lp.control(lpm$model)

  lpm$epsilon   <- sqrt(lpm$ctl$epsilon["epsint"])

  return(lpm)
}


###########################################################################################
#
# Calculate objfn values - Step 3
# Precompute but do not set objective function vectors. These are the same for every value of K,
# but in two phase model need to switch between two.
#
###########################################################################################

# Phase 1 - same for dea, fast, sdea - but if cook need to include M value
.dea_obj_phase1 <- function(lpm, vars.cook=FALSE, debug=1){

  obj.phase1                <- rep.int(0, lpm$col.n)
  obj.phase1[lpm$col.eff]   <- 1

  if(vars.cook){
    obj.phase1[lpm$col.exc] <- ifelse(lpm$orientation.in, 10^5, -10^5)
  }
  return(obj.phase1)
}

# Phase 2 - only used in two phase, values depend on if doing slack maximization or
#   secondary objective. Same for all models
.dea_obj_phase2 <- function(lpm, zT, vars.slack, debug=1){

  obj.phase2                <- rep.int(0, lpm$col.n)
  obj.phase2[1:lpm$nT]      <- ifelse( rep.int(lpm$orientation.in, lpm$nT), -zT, zT)

  # MUST be if, if slacks not set, slack vars may not exist in model - can not use ifelse
  if (vars.slack)                           # Maxamize Slacks  - set slacks to -1, 1
    obj.phase2[lpm$col.slack] <- ifelse(lpm$orientation.in, -1, 1)

  return(obj.phase2)
}


###########################################################################################
#
# Setup linear equations for RTS - Step 4
# Same for dea, sdea (with and without Cook), and fast
#
###########################################################################################

.dea_set_rts <- function(lpm, debug=1){

  # options.rts.l         <- c("vrs","drs", "crs", "irs")   # Set in global file
  rts.rhs               <- c(  1,   1,     0,     1)
  rts.typ               <- c( "=", "<=",   0,    ">=")

  # lookup rts in list of rts options and return index number
  rts.n = match(lpm$rts, options.rts.l)

  #  Scale constraints
  set.row(lpm$model, lpm$row.rts, rep.int(1, lpm$nT), c(1:lpm$nT))      # Fill in 1's for DMU's
  set.rhs(lpm$model,         rts.rhs[rts.n],  lpm$row.rts)              # 1 or 0
  set.constr.type(lpm$model, rts.typ[rts.n],  lpm$row.rts)              # =, >= or <=

  if (debug >= 3) lpDebugFile(lpm$model, "rts constraints added")

  return(lpm)

}

###########################################################################################
#
# Setup constant parts of technology set - Step 5
# Setup technology set parts of linear equation - X, Y, optionally slacks, constraints for variables
# These values are the same for each K
# Same for dea, sdea - if using slack vars need to setup slack vars
#
###########################################################################################

.dea_set_technology <- function(lpm, xT, yT, vars.slack=TRUE, debug=1){

  nx <- ncol(xT)
  ny <- ncol(yT)
  nv <- nx + ny

  # setup constraint for x's and y's - if slack vars =, otherwise
  constraint.type <- ifelse(vars.slack, "=", ">=")
  set.constr.type(lpm$model, rep.int(constraint.type, nv),    c(1:nv))

  # Setup constant x &  y's
  for (i in 1:nx) {                                         # Setup IN's (technology set)
    set.row(lpm$model, i,      -xT[, i],      c(1:lpm$nT))
  }

  for (r in 1:ny) {                                         # Setup OUT's (technology set)
    set.row(lpm$model, nx+r,    yT[, r],      c(1:lpm$nT))
  }

  # Setup slack varables - even if not using them - need to set so not free to be assigned
  if (vars.slack){
    col.slack <- lpm$col.slack[1] - 1
    for(i in 1:nv)
      set.mat(lpm$model,  i,  col.slack + i,  -1)              # Set Slack Diag to -1
  }

  if (debug >= 2) lpDebugFile(lpm$model, "technology part linear equation added")

  return(lpm)
}

###########################################################################################
#
# Precompute Values for inner loop for each K - Step 6
# To make inner loop for each K fast - all computations are done before loop
# The inner loop changes for each value of K being solved for.
#
###########################################################################################

# Compute efficiency values for column. Same for dea, sdea, fast
.dea_values_eff <- function(lpm, x, y, debug=1){

  nx <- ncol(x)
  ny <- ncol(y)
  nd <- nrow(x)

  test.x <- matrix(lpm$orientation.in, nrow=nd, ncol=nx)    # Ifelse produces output same
  test.y <- matrix(lpm$orientation.in, nrow=nd, ncol=ny)    # shape as logical variable

  # NOTE: setting column in model clears objfn for that colume. So we need to include objfn
  # value as row zero value at start of column

  objfn <- 1
  values.eff <- cbind( matrix(objfn, nrow=nd, ncol=1),     # Objfn value
                       ifelse(test.x,          x,    0),        ifelse(test.y,        0,    -y))

  return(values.eff)
}

# Return Excess Values - only used if computing cook sdea
.dea_values_exc <- function(lpm, x, y, debug=1){

  objfn <- ifelse(lpm$orientation.in, 10^5, -10^5)

  nx <- ncol(x)
  ny <- ncol(y)
  nd <- nrow(x)

  test.x <- matrix(lpm$orientation.in, nrow=nd, ncol=nx)  # Ifelse produces output same
  test.y <- matrix(lpm$orientation.in, nrow=nd, ncol=ny)  # shape as logical

  values.exc <- cbind( matrix(objfn, nrow=nd, ncol=1),     # Objfn value
                       ifelse(test.x,          0,    x),        ifelse(test.y,        y,    0))

  return(values.exc)
}

# Compute rhs values for column. Same for dea, sdea, fast - but if doing cook different
.dea_values_rhs <- function(lpm, x, y, cook=FALSE, debug=1){

  nx <- ncol(x)
  ny <- ncol(y)
  nd <- nrow(x)

  test.x <- matrix(lpm$orientation.in, nrow=nd, ncol=nx)    # Ifelse produces output same
  test.y <- matrix(lpm$orientation.in, nrow=nd, ncol=ny)    #  shape as logical

  values.rhs <- cbind(
    ifelse(cook | !test.x, -x,    0),        ifelse(cook | test.y,     y,     0))

  return(values.rhs)
}

# Compute super values - if doing sdea, fast(sdea) it forces DMU K lambda to be zero
# If not doing super, is entire array of zeros, Is complex because of index.K, index.T
# Possible options
# -K is not in technology set (index.T) - so no lambda needs to be zero
# -K is in technology set, need to figure out what col it is in
#
# ToDo: recode with no loop? lapply?
.dea_values_super <- function(lpm, nd, index.K, index.T, super, debug=1){

  values.super <- matrix(0, nrow=nd, ncol=lpm$col.n)
  match.K <- rep.int(0, nd)

  match.K[index.K] <-  match(index.K, index.T, nomatch=0)
  for(i in index.K){
    if (match.K[i] != 0)
      values.super[i, match.K[i]] <- 1 && super
  }

  return(values.super)
}

###########################################################################################
#
# Solve for Each K - Step 7
# Same for dea, sdea or fast. Same for phase 1 and phase 2
#
###########################################################################################

.dea_solveLPM <- function(lpm, k, phase, lpm.dual.n, debug=1){

  status <- solve.lpExtPtr(lpm$model)

  # If status is numerical instablity (5), set bias, try again, update vars - may work
  # Status 5 = "numerical failure encountered
  if (status == 5){
    set.basis(lpm$model, default = TRUE)
    status <- solve.lpExtPtr(lpm$model)
    if (debug >= 1) cat("Solver Phase ", phase, "First Try Status: ", status, "DMU k=", k,
                        "Retrying...\n")
  }

  # Need to rep status to make as wide as results vector
  # If status != 0, make results NA
  tmp.var <- ifelse( rep.int(status == 0, lpm$col.n),
                     get.variables(lpm$model)[1:lpm$col.n],
                     NA)

  # non-zero force excess to 0 not NA
  tmp.exc <- ifelse( status==0,
                     tmp.var[lpm$col.exc],
                     0)

  # for status 2 or 3 - make Eff +Inf or -Inf
  tmp.eff <- ifelse( status==0,
                     tmp.var[lpm$col.eff],
                     ifelse( status == 2 || status == 3,
                             ifelse(lpm$orientation.in, Inf, -Inf),
                             NA))

  tmp.dual <- ifelse( rep.int(status == 0, lpm.dual.n),
                      ifelse( rep.int(lpm$orientation.in, lpm.dual.n),
                              get.dual.solution(lpm$model), -get.dual.solution(lpm$model)),
                      NA)

  if (status == 2 || status == 3){
    msg <- lp_solve_error_msg(status)
    if (debug >= 2) cat("Solver Phase ", phase, "Status: ", status, "DMU k=", k,
                        " not in technology set: ", msg, "\n")
  } else if (status != 0){
    msg <- lp_solve_error_msg(status)
    if (debug >= 1) cat("Solver Phase ", phase, "Status: ", status, "DMU k=", k,
                        " encountered solver failure: ", msg, "\n")
  }

  tmp.var[lpm$col.eff] <- tmp.eff

  if(lpm$vars.cook){
    tmp.var[lpm$col.exc] <- tmp.exc
  }

  return (list(status=status, var=tmp.var, dual=tmp.dual))
}



