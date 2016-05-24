
######## input minimization with rescaling ######
dea.input.rescaling <- function(XREF, YREF, X, Y, RTS, rescaling=1) {
  XREF.rescaled = XREF * rescaling
  result = dea.input(XREF=XREF.rescaled, YREF=YREF, X=X, Y=Y, RTS=RTS)
  return(result$thetaOpt)
}

dea.output.rescaling <- function(XREF, YREF, X, Y, RTS, rescaling=1) {
  YREF.rescaled = YREF / rescaling
  result = dea.output(XREF=XREF, YREF=YREF.rescaled, X=X, Y=Y, RTS=RTS)
  return(result$thetaOpt)
}

dea.costmin.rescaling <- function(XREF, YREF, W, X, Y, RTS, rescaling=1) {
  XREF.rescaled = XREF * rescaling
  result = dea.costmin(XREF=XREF.rescaled, YREF=YREF, W=W, X=X, Y=Y, RTS=RTS)
  return( result$gammaOpt )
}


dea <- function(XREF, YREF, X, Y, W=NULL, model, RTS="variable") {
  if (missing(model) || ! model %in% c("input", "output", "costmin") ) {
    stop("Parameter 'model' has to be either 'input', 'output' or 'costmin'.")
  }
  if (missing(XREF)) stop("XREF is missing.")
  if (missing(YREF)) stop("YREF is missing.")
  if (missing(X))    stop("X is missing.")
  if (missing(Y))    stop("Y is missing.")
  if (model == "input") {
    return( dea.input(XREF=XREF, YREF=YREF, X=X, Y=Y, RTS=RTS) )
  } else if (model == "output") {
    return( dea.output(XREF=XREF, YREF=YREF, X=X, Y=Y, RTS=RTS) )
  } else { 
    ## costmin
    return( dea.costmin(XREF=XREF, YREF=YREF, X=X, Y=Y, W=W, RTS=RTS) )
  }
  
}

######## cost minimization with LP using GLPK ########
### XREF - N x Dinput matrix for ref. inputs (possibly multiplied by prices), N is the number of firms
### YREF - N x Doutput matrix for ref. outputs
### X - M x Dinput matrix of inputs (only for calculating cost efficiency)
### Y - M x Doutput matrix of outputs
### RTS - returns to scale: "variable" (default), "constant", "non-increasing"
######## output #######
### out$thetaOpt = thetaOpt
### out$lambda   = lambda
### out$feasible = feasible
### out$lambda_sum = rowSums(lambda)
dea.input <- function(XREF, YREF, X, Y, RTS="variable") {
  if ( ! is.matrix(XREF)) { XREF = as.matrix(XREF) }
  if ( ! is.matrix(YREF)) { YREF = as.matrix(YREF) }
  if ( ! is.matrix(X)) { X = as.matrix(X) }
  if ( ! is.matrix(Y)) { Y = as.matrix(Y) }

  if ( ! is.numeric(XREF)) { stop("XREF has to be numeric matrix or data.frame.") }
  if ( ! is.numeric(YREF)) { stop("YREF has to be numeric matrix or data.frame.") }
  if ( ! is.numeric(X)) { stop("X has to be numeric matrix or data.frame.") }
  if ( ! is.numeric(Y)) { stop("Y has to be numeric matrix or data.frame.") }
  
  if (nrow(XREF) != nrow(YREF)) { stop( sprintf("Number of rows in XREF (%d) does not equal the number of rows in YREF (%d)", nrow(XREF), nrow(YREF)) ) }
  if (nrow(X)    != nrow(Y)) { stop( sprintf("Number of rows in X (%d) does not equal the number of rows in Y (%d)", nrow(X), nrow(Y)) ) }
  if (ncol(XREF) != ncol(X)) { stop( sprintf("Number of columns in XREF (%d) does not equal the number of columns in X (%d)", ncol(XREF), ncol(X)) ) }
  if (ncol(YREF) != ncol(Y)) { stop( sprintf("Number of columns in YREF (%d) does not equal the number of columns in Y (%d)", ncol(YREF), ncol(Y)) ) }
  
  ## N is the number of firms in reference
  N = nrow(XREF)
  ## M is the number of firms for which we calculate the efficiency
  M = nrow(X)
  Dinput  = ncol(XREF)
  Doutput = ncol(YREF)
  
  ## adding 1 constraint if RTS is "variable" or "non-increasing"
  Drts = RTS %in% c("variable", "non-increasing")
  
  # variables = c(lambdas(N), theta(1) )
  
  # objective function:
  obj = c(rep(0, N), 1)
  
  # constraints in GLPK by default set variables to [0, Inf), see bounds in GLPK
  # constraint matrix C, RHS, constraint type:
  # we have a constriant for each output and for each input.
  C  = matrix(0.0, Doutput+Dinput+Drts, N+1 )
  b  = rep(0.0,    Doutput+Dinput+Drts)
  cd = rep(">=",   Doutput+Dinput+Drts)
  
  # constraints on outputs:
  C[1:Doutput, 1:N] = t(YREF)
  # b[1:Doutput] will be changed for each firm in the for loop
  
  # constraints on inputs: -X'*lambda + theta*x0 >= 0
  C[(Doutput+1):(Doutput+Dinput), 1:N ]  = -t(XREF)  ## -X' * lambda
  #C[(Doutput+1):(Doutput+Dinput), (N+1)] will be changed for each firm in the loop
  
  # positivity constrant on lambdas (unnecessary, they are bounded below from 0 in GLPK):
  
  # variable returns to scale constraint:
  if (RTS=="variable") {
    ## adding constraint sum(lambda) == 1
    i = Doutput + Dinput + 1
    C[i,] = c(rep(1.0,N), 0.0)
    b[i]  = 1.0
    cd[i] = "=="
  } else if (RTS=="non-increasing") {
    ## adding constraint sum(lambda) < 1
    i = Doutput + Dinput + 1
    C[i,] = c(rep(1.0,N), 0.0)
    b[i]  = 1.0
    cd[i] = "<"
  } else {
    ## constant returns to scale, no constraint needed.
  }
  
  ## using multi optimization to solve M linear programs
  outlp = multi_glpk_solve_LP(obj=obj, mat=C, dir=cd, rhs=b,
                              mrhs_i   = 1:Doutput,
                              mrhs_val = t(Y),
                              mmat_i   = cbind( (Doutput+1):(Doutput+Dinput), N+1 ),
                              mmat_val = t(X))
  # output:
  out = list()
  out$feasible = outlp$status==0
  out$thetaOpt = outlp$solution[N+1,]
  out$lambda   = t(outlp$solution[1:N, ])
  # necessary for VRS case, all rows should be equal to 1:
  out$lambda_sum = rowSums(out$lambda)
  
  return(out)
}


######## output DEA with GLPK ########
### XREF - N x Dinput matrix for ref. inputs (possibly multiplied by prices), N is the number of firms
### YREF - N x Doutput matrix for ref. outputs
### X - M x Dinput matrix of inputs (only for calculating cost efficiency)
### Y - M x Doutput matrix of outputs
### RTS - returns to scale: "variable" (default), "constant", "non-increasing"
######## output #######
### out$thetaOpt = thetaOpt
### out$lambda   = lambda
### out$feasible = feasible
### out$lambda_sum = rowSums(lambda)
dea.output <- function(XREF, YREF, X, Y, RTS="variable") {
  if ( ! is.matrix(XREF)) { XREF = as.matrix(XREF) }
  if ( ! is.matrix(YREF)) { YREF = as.matrix(YREF) }
  if ( ! is.matrix(X)) { X = as.matrix(X) }
  if ( ! is.matrix(Y)) { Y = as.matrix(Y) }
  
  if ( ! is.numeric(XREF)) { stop("XREF has to be numeric matrix or data.frame.") }
  if ( ! is.numeric(YREF)) { stop("YREF has to be numeric matrix or data.frame.") }
  if ( ! is.numeric(X)) { stop("X has to be numeric matrix or data.frame.") }
  if ( ! is.numeric(Y)) { stop("Y has to be numeric matrix or data.frame.") }
  
  if (nrow(XREF) != nrow(YREF)) { stop( sprintf("Number of rows in XREF (%d) does not equal the number of rows in YREF (%d)", nrow(XREF), nrow(YREF)) ) }
  if (nrow(X)    != nrow(Y)) { stop( sprintf("Number of rows in X (%d) does not equal the number of rows in Y (%d)", nrow(X), nrow(Y)) ) }
  if (ncol(XREF) != ncol(X)) { stop( sprintf("Number of columns in XREF (%d) does not equal the number of columns in X (%d)", ncol(XREF), ncol(X)) ) }
  if (ncol(YREF) != ncol(Y)) { stop( sprintf("Number of columns in YREF (%d) does not equal the number of columns in Y (%d)", ncol(YREF), ncol(Y)) ) }
  
  ## N is the number of firms in reference
  N = nrow(XREF)
  ## M is the number of firms for which we calculate the efficiency
  M = nrow(X)
  Dinput  = ncol(XREF)
  Doutput = ncol(YREF)
  
  ## adding 1 constraint if RTS is "variable" or "non-increasing"
  Drts = RTS %in% c("variable", "non-increasing")
  
  # variables = c(lambdas(N), theta(1) )
  
  # objective function:
  obj = c(rep(0, N), 1)
  
  # constraints in Rglpk by default set variables to [0, Inf)
  # constraint matrix C, RHS, constraint type:
  # we have a constriant for each output and for each input.
  C  = matrix(0.0, Doutput+Dinput+Drts, N+1 )
  b  = rep(0.0,    Doutput+Dinput+Drts)
  cd = rep(">=",   Doutput+Dinput+Drts)
  
  # constraints on outputs: Y'*lambda - theta*y >= 0
  C[1:Doutput, 1:N] = t(YREF)
  #C[1:Doutput, (N+1)] will be changed for each firm in the loop
  
  # constraints on inputs: -X'*lambda >= -x0
  C[(Doutput+1):(Doutput+Dinput), 1:N ]  = -t(XREF)  ## -X' * lambda
  # b[(Doutput+1):(Doutput+Dinput)] will be changed for each firm in the for loop
  
  # positivity constrant on lambdas (unnecessary, they are bounded below from 0 in Rglpk):
  
  # variable returns to scale constraint:
  if (RTS=="variable") {
    ## adding constraint sum(lambda) == 1
    i = Doutput + Dinput + 1
    C[i,] = c(rep(1.0,N), 0.0)
    b[i]  = 1.0
    cd[i] = "=="
  } else if (RTS=="non-increasing") {
    ## adding constraint sum(lambda) < 1
    i = Doutput + Dinput + 1
    C[i,] = c(rep(1.0,N), 0.0)
    b[i]  = 1.0
    cd[i] = "<"
  } else {
    ## constant returns to scale, no constraint needed.
  }

  outlp = multi_glpk_solve_LP(obj=obj, mat=C, dir=cd, rhs=b,
                              mrhs_i   = (Doutput+1):(Doutput+Dinput),
                              mrhs_val = t(-X),
                              mmat_i   = cbind(1:Doutput, N+1),
                              mmat_val = t(-Y),
                              max=TRUE )
  # output:
  out = list()
  out$feasible = outlp$status==0
  out$thetaOpt = 1 / outlp$solution[N+1,]
  out$lambda   = t(outlp$solution[1:N, ])
  # necessary for VRS case, all rows should be equal to 1:
  out$lambda_sum = rowSums(out$lambda)
  
  return(out)
}


######## cost minimization with LP using GLPK ########
### XREF - N x Dinput matrix for ref. inputs (possibly multiplied by prices), N is the number of firms
### YREF - N x Doutput matrix for ref. outputs
### X    - M x Dinput matrix of inputs (only for calculating cost efficiency)
### Y    - M x Doutput matrix of outputs
### W    - M x Dinput matrix of input prices (ones for Tone version), M is number of investigated firms
### RTS  - returns to scale, "variable" (default), "constant" and "non-increasing"
dea.costmin <- function(XREF, YREF, X, Y, W, RTS="variable") {
  if ( ! is.matrix(XREF)) { XREF = as.matrix(XREF) }
  if ( ! is.matrix(YREF)) { YREF = as.matrix(YREF) }
  if ( ! is.matrix(X)) { X = as.matrix(X) }
  if ( ! is.matrix(Y)) { Y = as.matrix(Y) }
  if ( ! is.matrix(W)) { W = as.matrix(W) }
  
  if ( ! is.numeric(XREF)) { stop("XREF has to be numeric matrix or data.frame.") }
  if ( ! is.numeric(YREF)) { stop("YREF has to be numeric matrix or data.frame.") }
  if ( ! is.numeric(X)) { stop("X has to be numeric matrix or data.frame.") }
  if ( ! is.numeric(Y)) { stop("Y has to be numeric matrix or data.frame.") }
  if ( ! is.numeric(W)) { stop("W has to be numeric matrix or data.frame.") }
  
  if (nrow(XREF) != nrow(YREF)) { stop( sprintf("Number of rows in XREF (%d) does not equal the number of rows in YREF (%d)", nrow(XREF), nrow(YREF)) ) }
  if (nrow(X)    != nrow(Y)) { stop( sprintf("Number of rows in X (%d) does not equal the number of rows in Y (%d)", nrow(X), nrow(Y)) ) }
  if (ncol(XREF) != ncol(X)) { stop( sprintf("Number of columns in XREF (%d) does not equal the number of columns in X (%d)", ncol(XREF), ncol(X)) ) }
  if (ncol(YREF) != ncol(Y)) { stop( sprintf("Number of columns in YREF (%d) does not equal the number of columns in Y (%d)", ncol(YREF), ncol(Y)) ) }
  if (nrow(X)    != nrow(W)) { stop( sprintf("Number of rows in X (%d) does not equal the number of rows in W (%d)", nrow(X), nrow(W)) ) }
  if (ncol(X)    != ncol(W)) { stop( sprintf("Number of columns in X (%d) does not equal the number of columns in W (%d)", ncol(X), ncol(W)) ) }
  
  N = nrow(XREF)
  M = nrow(W)
  Dinput  = ncol(XREF)
  Doutput = ncol(YREF)

  ## adding 1 constraint if RTS is "variable" or "non-increasing"
  Drts = RTS %in% c("variable", "non-increasing")
  
  # variables = c(lambdas(N), x(Dinput) )
  
  # objective function:
  obj = rep.int(0.0, N + Dinput)
  
  # constraints in GLPK by default set variables to [0, Inf), see bounds in GLPK
  # constraint matrix, RHS, constraint type:
  C  = matrix(0.0, Doutput+Dinput+Drts, N+Dinput )
  b  = rep(0.0,    Doutput+Dinput+Drts)
  cd = rep(">=",   Doutput+Dinput+Drts)
  
  # constraints on outputs:
  C[1:Doutput, 1:N] = t(YREF)
  #b[1:Doutput] = Y
  
  # constraints on inputs:
  C[(Doutput+1):(Doutput+Dinput), 1:N ]             = -t(XREF)
  C[(Doutput+1):(Doutput+Dinput), (N+1):(N+Dinput)] = diag(Dinput)
  
  # constrant on lambdas (unnecessary, they are bounded below from 0 in Rglpk):
  #n = Doutput+Dinput
  #C[(n+1):(n+N), 1:N] = diag(N)
  
  # variable returns to scale constraint:
  if (RTS=="variable") {
    i = Doutput + Dinput + 1
    C[i,] = rep.int(c(1.0, 0.0), c(N, Dinput))
    b[i]  = 1.0
    cd[i] = "=="
  } else if (RTS=="non-increasing") {
    i = Doutput + Dinput + 1
    C[i,] = rep.int(c(1.0, 0.0), c(N, Dinput))
    b[i]  = 1.0
    cd[i] = "<"
  }
  
  outlp = multi_glpk_solve_LP(obj=obj, mat=C, dir=cd, rhs=b,
                              mobj_i   = (N+1):(N+Dinput),
                              mobj_val = t(W),
                              mrhs_i   = 1:Doutput,
                              mrhs_val = t(Y))
  # output
  out = list()
  out$XOpt     = t(outlp$solution[(N+1):(N+Dinput), , drop=F])
  out$gammaOpt = unname( rowSums(out$XOpt*W) / rowSums(X*W) )
  out$feasible = outlp$status==0
  out$lambda   = t(outlp$solution[1:N, ])
  
  # necessary for VRS case, all rows should be equal to 1:
  out$lambda_sum = rowSums(out$lambda)
  
  return(out)
}
