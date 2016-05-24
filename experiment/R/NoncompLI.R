NoncompLI <- function(formulae, Z, D, data = parent.frame(), n.draws = 5000,
                      param = TRUE, in.sample = FALSE, model.c = "probit",
                      model.o = "probit", model.r = "probit", 
                      tune.c = 0.01, tune.o = 0.01, tune.r = 0.01,
                      tune.v = 0.01, p.mean.c = 0, p.mean.o = 0,
                      p.mean.r = 0, p.prec.c = 0.001,
                      p.prec.o = 0.001, p.prec.r = 0.001,
                      p.df.o = 10, p.scale.o = 1, p.shape.o = 1,
                      mda.probit = TRUE, coef.start.c = 0,
                      coef.start.o = 0, tau.start.o = NULL,
                      coef.start.r = 0, var.start.o = 1,
                      burnin = 0, thin = 0, verbose = TRUE) {  

  ## getting the data
  call <- match.call()

  ## model types
  if (!(model.c %in% c("logit", "probit"))) 
    stop("no such model is supported for the compliance model.")
  
  if (!(model.o %in% c("logit", "probit", "oprobit", "gaussian",
                       "negbin", "twopart"))) 
    stop("no such model is supported for the outcome model.")
  
  if (!(model.r %in% c("logit", "probit"))) 
    stop("no such model is supported for the response model.")    

  ## outcome model
  mf <- model.frame(formulae[[1]], data=data, na.action='na.pass')
  Xo <- model.matrix(formulae[[1]], data=mf)
  if (sum(is.na(Xo)) > 0)
    stop("missing values not allowed in covariates")
  if (model.o %in% c("gaussian", "twopart")) {
    Y <- model.response(mf)
    if (model.o == "twopart") {
      Y1 <- Y
      Y[!is.na(Y)] <- (Y[!is.na(Y)] > 0)*1
    }
  } else if (model.o == "oprobit")
    Y <- as.integer(factor(model.response(mf)))
  else
    Y <- as.integer(model.response(mf))

  ## compliance model
  Xc <- model.matrix(formulae[[2]], data=data)

  ## response model
  if (any(is.na(Y)))
    Xr <- model.matrix(formulae[[3]], data=data)
  else
    Xr <- model.matrix(~ 1, data = data)
  
  N <- length(Y)
  Z <- eval(call$Z, envir = data)
  D <- eval(call$D, envir = data)
  if (sum(is.na(Z)) > 0)
    stop("missing values not allowed in the encouragement variable")

  res <- list(call = call, Y = Y, Xo = Xo, Xc = Xc, Xr = Xr,
              D = D, Z = Z, n.draws = n.draws)
  if (model.o == "twopart") {
    res$Y1 <- Y1
    nsamp1 <- sum(Y1[!is.na(Y1)] > 0)
    Y1[is.na(Y1)] <- 0
  }
  
  ## Starting values for missing D
  RD <- (!is.na(D))*1
  NRD <- is.na(D)
  if (sum(NRD) > 0)
    D[NRD] <- Z[NRD]

  ## Random starting values for missing Y using Bernoulli(0.5) for
  ## binary and ordinal
  R <- (!is.na(Y))*1
  NR <- is.na(Y)
  Ymiss <- sum(NR)
  if (Ymiss > 0)
    if (model.o == "gaussian")
      Y[NR] <- rnorm(Ymiss)
    else
      Y[NR] <- (runif(Ymiss) > 0.5)*1
  if (model.o == "oprobit") {
    ncat <- max(Y, na.rm = TRUE) + 1
    if (is.null(tau.start.o))
      tau.start.o <- seq(from = 0, length = ncat-1)/10
    if (length(tau.start.o) != (ncat-1))
      stop("incorrect length for tau.start.o")
    if (!identical(sort(tau.start.o), tau.start.o))
      stop("incorrect input for tau.start.o")
    if (length(unique(tau.start.o)) != (ncat-1))
      stop("incorrect input for tau.start.o")
    tau.start.o <- c(tau.start.o, tau.start.o[ncat-1]+1000)
  }

  ## Compliance status: 0 = noncomplier, 1 = complier
  C <- rep(NA, N)
  C[Z == 1 &  D == 0] <- 0 # never-takers

  ## Always-takers: 0 = never-takers/compliers, 1 = always-takers
  if (sum(Z == 0 & D == 1)>0) { # some always-takers
    AT <- TRUE
    A <- rep(NA, N)
    A[D == 0] <- 0            # never-takers or compliers
    A[Z == 0 & D == 1] <- 1  
    if (model.c == "logit")
      C[Z == 0 & D == 1] <- 2 # always-takers
    else
      C[Z == 0 & D == 1] <- 0 # always-takers
  } else { # no always-takers
    A <- rep(0, N)
    AT <- FALSE
    C[Z == 1 & D == 1] <- 1   # compliers
  }
  res$R <- R
  res$A <- A
  res$C <- C
  
  ## Random starting values for missing compliance status
  if (AT) {
    A[is.na(A)] <- (runif(sum(is.na(A))) > 0.5)*1
    if (model.c == "logit")
      C[A == 1] <- 2
    else
      C[A == 1] <- 0
  }
  C[is.na(C)] <- (runif(sum(is.na(C))) > 0.5)*1
  
  ## Completing the outcome model matrix 
  ## The default category is never-takers
  X <- Xo
  X1 <- Xr
  if (AT) { # when some always-takers exist
  ## Xo = [c1 c0 a X] where c1 for compliers with encouragement
  ##                        c0 for compliers without encouragement
  ##                        a for always-takers with/without encouragement  
    Xo <- cbind(0, 0, 0, X)
    Xr <- cbind(0, 0, 0, X1)
    Xo[A == 1, 3] <- 1
    Xr[A == 1, 3] <- 1
    colnames(Xo) <- c("Complier1", "Complier0", "AlwaysTaker",
                      colnames(X))
    colnames(Xr) <- c("Complier1", "Complier0", "AlwaysTaker",
                      colnames(X1))
  } else { # when always-takers do not exist
  ## Xo = [c1 c0 X] where c1 for compliers with encouragement
  ##                      c0 for compliers without encouragement
    Xo <- cbind(0, 0, X)
    Xr <- cbind(0, 0, X1)
    colnames(Xo) <- c("Complier1", "Complier0", colnames(X))
    colnames(Xr) <- c("Complier1", "Complier0", colnames(X1))
  }
  Xo[C == 1 & Z == 1, 1] <- 1
  Xo[C == 1 & Z == 0, 2] <- 1
  Xr[C == 1 & Z == 1, 1] <- 1
  Xr[C == 1 & Z == 0, 2] <- 1
  
  ## dimensions
  ncovC <- ncol(Xc)
  ncovO <- ncol(Xo)
  ncovR <- ncol(Xr)
  if (model.o == "oprobit")
    if (AT)
      nqoi <- 2 + (ncat-1)*6
    else
      nqoi <- 2 + (ncat-1)*5
  else
    if (AT)
      nqoi <- 8
    else
      nqoi <- 7

  ## checking starting values and prior
  if ((model.c == "logit") & AT) {
    if(length(p.mean.c) != ncovC*2)
      if (length(p.mean.c) == 1)
        p.mean.c <- rep(p.mean.c, ncovC*2)
      else
        stop(paste("the length of p.mean.c should be", ncovC*2))    
    if(length(coef.start.c) != ncovC*2)
      if (length(coef.start.c) == 1)
        coef.start.c <- rep(coef.start.c, ncovC*2)
      else
        stop(paste("the length of coef.start.c should be", ncovC*2))        
  } else {
    if(length(p.mean.c) != ncovC)
      if (length(p.mean.c) == 1)
        p.mean.c <- rep(p.mean.c, ncovC)
      else
        stop(paste("the length of p.mean.c should be", ncovC))        
    if(length(coef.start.c) != ncovC)
      if (length(coef.start.c) == 1)
        coef.start.c <- rep(coef.start.c, ncovC)
      else
        stop(paste("the length of coef.start.c should be", ncovC))    
  }

  if(length(coef.start.o) != ncovO)
    if (length(coef.start.o) == 1)
      coef.start.o <- rep(coef.start.o, ncovO)
    else
      stop(paste("the length of coef.start.o should be", ncovO))      
  if(length(p.mean.o) != ncovO)
    if (length(p.mean.o) == 1)
      p.mean.o <- rep(p.mean.o, ncovO)
    else
      stop(paste("the length of p.mean.o should be", ncovO))    

  if(length(coef.start.r) != ncovR)
    if (length(coef.start.r) == 1)
      coef.start.r <- rep(coef.start.r, ncovR)
    else
      stop(paste("the length of coef.start.r should be", ncovR))    
  if(length(p.mean.r) != ncovR)
    if (length(p.mean.r) == 1)
      p.mean.r <- rep(p.mean.r, ncovR)
    else
      stop(paste("the length of p.mean.r should be", ncovR))    

  if(is.matrix(p.prec.c)) {
    if (dim(p.prec.c) != rep(ncovC*2, 2))
        stop(paste("the dimension of p.prec.c should be",
                   rep(ncovC*2, 2)))    
  } else if (length(p.prec.c) == 1){
    if (model.c == "logit" & AT)
      p.prec.c <- diag(p.prec.c, ncovC*2)
    else
      p.prec.c <- diag(p.prec.c, ncovC)
  } else {
    stop("Incorrect input for p.prec.c")
  }

  if(is.matrix(p.prec.o)) {
    if (dim(p.prec.o) != rep(ncovO, 2))
      stop(paste("the dimension of p.prec.o should be",
                 rep(ncovO, 2)))    
  } else if (length(p.prec.o) == 1){
    p.prec.o <- diag(p.prec.o, ncovO)
  } else {
    stop("Incorrect input for p.prec.o")
  }

  if(is.matrix(p.prec.r)) {
    if (dim(p.prec.r) != rep(ncovR, 2))
      stop(paste("the dimension of p.prec.r should be",
                 rep(ncovR, 2)))    
  } else if (length(p.prec.r) == 1){
    p.prec.r <- diag(p.prec.r, ncovR)
  } else {
    stop("Incorrect input for p.prec.r")
  }
  
  ## proposal variance for logits
  if (model.c == "logit")
    if (AT) {
      if (length(tune.c) != ncovC*2)
        if (length(tune.c) == 1)
          tune.c <- rep(tune.c, ncovC*2)
        else
          stop(paste("the length of tune.c should be", ncovC*2))
    } else {
      if (length(tune.c) != ncovC)
        if (length(tune.c) == 1)
          tune.c <- rep(tune.c, ncovC)
        else
          stop(paste("the length of tune.c should be", ncovC))
    }
  if (model.o == "logit" || model.o == "negbin")
    if (length(tune.o) != ncovO)
      if (length(tune.o) == 1)
        tune.o <- rep(tune.o, ncovO)
      else
        stop(paste("the length of tune.o should be", ncovO))
  if (model.o == "oprobit")
    if (length(tune.o) != ncat-2)
      if (length(tune.o) == 1)
        tune.o <- rep(tune.o, ncat-2)
      else
        stop(paste("the length of tune.o should be", ncat-2))
  if (model.r == "logit")
    if (length(tune.r) != ncovR)
      if (length(tune.r) == 1)
        tune.r <- rep(tune.r, ncovR)
      else
        stop(paste("the length of tune.r should be", ncovR))
  
  ## checking thinnig and burnin intervals
  if (n.draws <= 0)
    stop("`n.draws' should be a positive integer.")
  if (burnin < 0 || burnin >= n.draws)
    stop("`burnin' should be a non-negative integer less than `n.draws'.")
  if (thin < 0 || thin >= n.draws)
    stop("`thin' should be a non-negative integer less than `n.draws'.")
  keep <- thin + 1

  ## calling C function
  if (model.o == "probit" || model.o == "logit")
    out <- .C("LIbinary",
              as.integer(Y), as.integer(R), as.integer(Z),
              as.integer(D), as.integer(RD), as.integer(C), as.integer(A),
              as.integer(Ymiss), as.integer(AT), as.integer(in.sample), 
              as.double(Xc), as.double(Xo), as.double(Xr),
              as.double(coef.start.c), as.double(coef.start.c),
              as.double(coef.start.o), as.double(coef.start.r),
              as.integer(N), as.integer(n.draws),
              as.integer(ncovC), as.integer(ncovO), as.integer(ncovR),
              as.double(p.mean.c), as.double(p.mean.o),
              as.double(p.mean.r),
              as.double(p.prec.c), as.double(p.prec.o),
              as.double(p.prec.r),
              as.double(tune.c), as.double(tune.o), as.double(tune.r),
              as.integer(model.c == "logit"),
              as.integer(model.o == "logit"),
              as.integer(model.r == "logit"),
              as.integer(param), as.integer(mda.probit), as.integer(burnin),
              as.integer(keep), as.integer(verbose),
              coefC = double(ncovC*(ceiling((n.draws-burnin)/keep))),
              coefA = double(ncovC*(ceiling((n.draws-burnin)/keep))),
              coefO = double(ncovO*(ceiling((n.draws-burnin)/keep))),
              coefR = double(ncovR*(ceiling((n.draws-burnin)/keep))),
              QoI = double(nqoi*(ceiling((n.draws-burnin)/keep))),
              PACKAGE = "experiment")
  else if (model.o == "oprobit")
    out <- .C("LIordinal",
              as.integer(Y), as.integer(R), as.integer(Z),
              as.integer(D), as.integer(RD), as.integer(C), as.integer(A),
              as.integer(Ymiss), as.integer(AT), as.integer(in.sample), 
              as.double(Xc), as.double(Xo), as.double(Xr),
              as.double(coef.start.c), as.double(coef.start.c),
              as.double(coef.start.o), as.double(tau.start.o),
              as.double(coef.start.r), 
              as.integer(N), as.integer(n.draws), as.integer(ncat),
              as.integer(ncovC), as.integer(ncovO), as.integer(ncovR),
              as.double(p.mean.c), as.double(p.mean.o),
              as.double(p.mean.r),
              as.double(p.prec.c), as.double(p.prec.o),
              as.double(p.prec.r),
              as.double(tune.c), as.double(tune.o), as.double(tune.r),
              as.integer(model.c == "logit"),
              as.integer(model.r == "logit"),
              as.integer(param), as.integer(mda.probit), as.integer(burnin),
              as.integer(keep), as.integer(verbose),
              coefC = double(ncovC*(ceiling((n.draws-burnin)/keep))),
              coefA = double(ncovC*(ceiling((n.draws-burnin)/keep))),
              coefO = double(ncovO*(ceiling((n.draws-burnin)/keep))),
              coefR = double(ncovR*(ceiling((n.draws-burnin)/keep))),
              tauO = double((ncat-1)*(ceiling((n.draws-burnin)/keep))),
              QoI = double(nqoi*(ceiling((n.draws-burnin)/keep))),
              PACKAGE = "experiment")
  else if (model.o == "gaussian")
    out <- .C("LIgaussian",
              as.double(Y), as.integer(R), as.integer(Z),
              as.integer(D), as.integer(RD), as.integer(C), as.integer(A),
              as.integer(Ymiss), as.integer(AT), as.integer(in.sample), 
              as.double(Xc), as.double(Xo), as.double(Xr),
              as.double(coef.start.c), as.double(coef.start.c),
              as.double(coef.start.o), as.double(var.start.o),
              as.double(coef.start.r),
              as.integer(N), as.integer(n.draws),
              as.integer(ncovC), as.integer(ncovO), as.integer(ncovR),
              as.double(p.mean.c), as.double(p.mean.o),
              as.double(p.mean.r), 
              as.double(p.prec.c), as.double(p.prec.o),
              as.double(p.prec.r), as.integer(p.df.o),
              as.double(p.scale.o),
              as.double(tune.c), as.double(tune.r),
              as.integer(model.c == "logit"),
              as.integer(model.r == "logit"),
              as.integer(param), as.integer(mda.probit), as.integer(burnin),
              as.integer(keep), as.integer(verbose),
              coefC = double(ncovC*(ceiling((n.draws-burnin)/keep))),
              coefA = double(ncovC*(ceiling((n.draws-burnin)/keep))),
              coefO = double(ncovO*(ceiling((n.draws-burnin)/keep))),
              coefR = double(ncovR*(ceiling((n.draws-burnin)/keep))),
              var = double(ceiling((n.draws-burnin)/keep)),
              QoI = double(nqoi*(ceiling((n.draws-burnin)/keep))),
              PACKAGE = "experiment")
  else if (model.o == "negbin")
    out <- .C("LIcount",
              as.integer(Y), as.integer(R), as.integer(Z),
              as.integer(D), as.integer(RD), as.integer(C), as.integer(A),
              as.integer(Ymiss), as.integer(AT), as.integer(in.sample), 
              as.double(Xc), as.double(Xo), as.double(Xr),
              as.double(coef.start.c), as.double(coef.start.c),
              as.double(coef.start.o), as.double(var.start.o),
              as.double(coef.start.r),
              as.integer(N), as.integer(n.draws),
              as.integer(ncovC), as.integer(ncovO), as.integer(ncovR),
              as.double(p.mean.c), as.double(p.mean.o),
              as.double(p.mean.r), 
              as.double(p.prec.c), as.double(p.prec.o),
              as.double(p.prec.r), as.double(p.shape.o),
              as.double(p.scale.o), as.double(tune.c),
              as.double(tune.r), as.double(tune.o), as.double(tune.v),
              as.integer(model.c == "logit"),
              as.integer(model.r == "logit"),
              as.integer(param), as.integer(mda.probit), as.integer(burnin),
              as.integer(keep), as.integer(verbose),
              coefC = double(ncovC*(ceiling((n.draws-burnin)/keep))),
              coefA = double(ncovC*(ceiling((n.draws-burnin)/keep))),
              coefO = double(ncovO*(ceiling((n.draws-burnin)/keep))),
              coefR = double(ncovR*(ceiling((n.draws-burnin)/keep))),
              var = double(ceiling((n.draws-burnin)/keep)),
              QoI = double(nqoi*(ceiling((n.draws-burnin)/keep))),
              PACKAGE = "experiment")
  else if (model.o == "twopart")
    out <- .C("LItwopart",
              as.integer(Y), as.double(Y1), as.integer(R), as.integer(Z),
              as.integer(D), as.integer(RD), as.integer(C), as.integer(A),
              as.integer(Ymiss), as.integer(AT), as.integer(in.sample), 
              as.double(Xc), as.double(Xo), as.double(Xr),
              as.double(coef.start.c), as.double(coef.start.c),
              as.double(coef.start.o), as.double(coef.start.o),
              as.double(var.start.o), as.double(coef.start.r),
              as.integer(c(N, nsamp1)), as.integer(n.draws),
              as.integer(ncovC), as.integer(ncovO), as.integer(ncovR),
              as.double(p.mean.c), as.double(p.mean.o),
              as.double(p.mean.r), 
              as.double(p.prec.c), as.double(p.prec.o),
              as.double(p.prec.r), as.integer(p.df.o),
              as.double(p.scale.o),
              as.double(tune.c), as.double(tune.r),
              as.integer(model.c == "logit"),
              as.integer(model.r == "logit"),
              as.integer(param), as.integer(mda.probit), as.integer(burnin),
              as.integer(keep), as.integer(verbose),
              coefC = double(ncovC*(ceiling((n.draws-burnin)/keep))),
              coefA = double(ncovC*(ceiling((n.draws-burnin)/keep))),
              coefO = double(ncovO*(ceiling((n.draws-burnin)/keep))),
              coefO1 = double(ncovO*(ceiling((n.draws-burnin)/keep))),
              coefR = double(ncovR*(ceiling((n.draws-burnin)/keep))),
              var = double(ceiling((n.draws-burnin)/keep)),
              QoI = double(nqoi*(ceiling((n.draws-burnin)/keep))),
              PACKAGE = "experiment")
  
  if (param) {
    res$coefC <- matrix(out$coefC, byrow = TRUE, ncol = ncovC)
    colnames(res$coefC) <- colnames(Xc)
    if (AT) {
      res$coefA <- matrix(out$coefA, byrow = TRUE, ncol = ncovC)
      colnames(res$coefA) <- colnames(Xc)
    }
    res$coefO <- matrix(out$coefO, byrow = TRUE, ncol = ncovO)
    colnames(res$coefO) <- colnames(Xo)
    if (model.o == "twopart") {
      res$coefO1 <- matrix(out$coefO1, byrow = TRUE, ncol = ncovO)
      colnames(res$coefO1) <- colnames(Xo)
    }
    if (model.o == "oprobit")
      res$tau <- matrix(out$tauO, byrow = TRUE, ncol = ncat - 1)
    if (Ymiss > 0) {
      res$coefR <- matrix(out$coefR, byrow = TRUE, ncol = ncovR)
      colnames(res$coefR) <- colnames(Xr)
    }
    if (model.o %in% c("gaussian", "negbin", "twopart"))
      res$sig2 <- out$var
  }

  QoI <- matrix(out$QoI, byrow = TRUE, ncol = nqoi)
  if (model.o == "oprobit") {
    res$ITT <- QoI[,1:(ncat-1)]
    res$CACE <- QoI[,ncat:(2*(ncat-1))]
    res$Y1barC <- QoI[,(2*(ncat-1)+1):(3*(ncat-1))]
    res$Y0barC <- QoI[,(3*(ncat-1)+1):(4*(ncat-1))]
    res$YbarN <- QoI[,(4*(ncat-1)+1):(5*(ncat-1))]
    res$pC <- QoI[,(5*(ncat-1)+1)]
    res$pN <- QoI[,(5*(ncat-1)+2)]
    if (AT) 
      res$YbarA <- QoI[,(5*(ncat-1)+3):(6*(ncat-1)+2)]
  } else {
    res$ITT <- QoI[,1]
    res$CACE <- QoI[,2]
    res$pC <- QoI[,3]
    res$pN <- QoI[,4]
    res$Y1barC <- QoI[,5]
    res$Y0barC <- QoI[,6]
    res$YbarN <- QoI[,7]
    if (AT) 
      res$YbarA <- QoI[,8]
  }
  if (AT) 
    res$pA <- 1-res$pC-res$pN

  class(res) <- "NoncompLI"
  return(res)
}
