ecoNP <- function(formula, data = parent.frame(), N = NULL, supplement = NULL,
                  context = FALSE, mu0 = 0, tau0 = 2, nu0 = 4, S0 = 10,
                  alpha = NULL, a0 = 1, b0 = 0.1, parameter = FALSE,
                  grid = FALSE, n.draws = 5000, burnin = 0, thin = 0,
                  verbose = FALSE){ 

 ## contextual effects
  if (context)
    ndim <- 3
  else
    ndim <- 2

  ## checking inputs
  if (burnin >= n.draws)
    stop("n.draws should be larger than burnin")

  if (length(mu0)==1)
    mu0 <- rep(mu0, ndim)
  else if (length(mu0)!=ndim)
    stop("invalid inputs for mu0")
  if (is.matrix(S0)) {
    if (any(dim(S0)!=ndim))
      stop("invalid inputs for S0")
  }
  else
    S0 <- diag(S0, ndim)
  
  mf <- match.call()
  
  ## getting X, Y and N
  tt <- terms(formula)
  attr(tt, "intercept") <- 0
  if (is.matrix(eval.parent(mf$data)))
    data <- as.data.frame(data)
  X <- model.matrix(tt, data)
  Y <- model.response(model.frame(tt, data = data))
  N <- eval(mf$N, data)

  ## alpha
  if (is.null(alpha)) {
    alpha.update <- TRUE
    alpha <- 0
  }
  else
    alpha.update <- FALSE

  ## checking the data and calculating the bounds 
  tmp <- checkdata(X, Y, supplement, ndim)
  bdd <- ecoBD(formula, data=data)
  W1min <- bdd$Wmin[order(tmp$order.old)[1:nrow(tmp$d)],1,1]
  W1max <- bdd$Wmax[order(tmp$order.old)[1:nrow(tmp$d)],1,1]
 
  ## fitting the model
  n.store <- floor((n.draws-burnin)/(thin+1))
  unit.par <- unit.w <- tmp$n.samp+tmp$samp.X1+tmp$samp.X0
  n.par <- n.store * unit.par
  n.w <- n.store * unit.w
  unit.a <- 1

  if (context) 
    res <- .C("cDPecoX", as.double(tmp$d), as.integer(tmp$n.samp),
              as.integer(n.draws), as.integer(burnin), as.integer(thin+1),
              as.integer(verbose), as.integer(nu0), as.double(tau0),
              as.double(mu0), as.double(S0), as.double(alpha),
              as.integer(alpha.update), as.double(a0), as.double(b0),
              as.integer(tmp$survey.yes), as.integer(tmp$survey.samp),
              as.double(tmp$survey.data), as.integer(tmp$X1type),
              as.integer(tmp$samp.X1), as.double(tmp$X1.W1),
              as.integer(tmp$X0type), as.integer(tmp$samp.X0),
              as.double(tmp$X0.W2), 
              as.double(W1min), as.double(W1max), 
              as.integer(parameter), as.integer(grid),
              pdSMu0=double(n.par), pdSMu1=double(n.par),
              pdSMu2=double(n.par),	
              pdSSig00=double(n.par), pdSSig01=double(n.par),
              pdSSig02=double(n.par), pdSSig11=double(n.par),
              pdSSig12=double(n.par), pdSSig22=double(n.par), 
              pdSW1=double(n.w), pdSW2=double(n.w), 
              pdSa=double(n.store), pdSn=integer(n.store), PACKAGE="eco")
  else 
    res <- .C("cDPeco", as.double(tmp$d), as.integer(tmp$n.samp),
              as.integer(n.draws), as.integer(burnin), as.integer(thin+1),
              as.integer(verbose), as.integer(nu0), as.double(tau0),
              as.double(mu0), as.double(S0), as.double(alpha),
              as.integer(alpha.update), as.double(a0), as.double(b0),
              as.integer(tmp$survey.yes), as.integer(tmp$survey.samp),
              as.double(tmp$survey.data), as.integer(tmp$X1type),
              as.integer(tmp$samp.X1), as.double(tmp$X1.W1),
              as.integer(tmp$X0type), as.integer(tmp$samp.X0),
              as.double(tmp$X0.W2), 
              as.double(W1min), as.double(W1max), 
              as.integer(parameter), as.integer(grid),
              pdSMu0=double(n.par), pdSMu1=double(n.par),
              pdSSig00=double(n.par), pdSSig01=double(n.par),
              pdSSig11=double(n.par), pdSW1=double(n.w), pdSW2=double(n.w), 
              pdSa=double(n.store), pdSn=integer(n.store), PACKAGE="eco")
  
  ## output
  W1.post <- matrix(res$pdSW1, n.store, unit.w, byrow=TRUE)[,tmp$order.old]
  W2.post <- matrix(res$pdSW2, n.store, unit.w, byrow=TRUE)[,tmp$order.old]
  W <- array(rbind(W1.post, W2.post), c(n.store, 2, unit.w))
  colnames(W) <- c("W1", "W2")
  res.out <- list(call = mf, X = X, Y = Y, N = N, W = W,
                  Wmin = bdd$Wmin[,1,], Wmax = bdd$Wmax[,1,],
                  burin = burnin, thin = thin, nu0 = nu0, tau0 = tau0,
                  mu0 = mu0, a0 = a0, b0 = b0, S0 = S0)

  ## optional outputs
  if (parameter){
    if (context) {
      mu1.post <- matrix(res$pdSMu0, n.store, unit.par, byrow=TRUE)[,tmp$order.old]
      mu2.post <- matrix(res$pdSMu1, n.store, unit.par, byrow=TRUE)[,tmp$order.old]
      mu3.post <- matrix(res$pdSMu2, n.store, unit.par, byrow=TRUE)[,tmp$order.old]
      Sigma11.post <- matrix(res$pdSSig00, n.store, unit.par, byrow=TRUE)[,tmp$order.old]
      Sigma12.post <- matrix(res$pdSSig01, n.store, unit.par, byrow=TRUE)[,tmp$order.old]
      Sigma13.post <- matrix(res$pdSSig02, n.store, unit.par, byrow=TRUE)[,tmp$order.old]
      Sigma23.post <- matrix(res$pdSSig12, n.store, unit.par, byrow=TRUE)[,tmp$order.old]
      Sigma22.post <- matrix(res$pdSSig11, n.store, unit.par, byrow=TRUE)[,tmp$order.old]
      Sigma33.post <- matrix(res$pdSSig22, n.store, unit.par, byrow=TRUE)[,tmp$order.old]
      res.out$mu <- array(rbind(mu1.post, mu2.post, mu3.post),
                          dim=c(n.store, 3, unit.par),
                          dimnames=list(1:n.store, c("mu1", "mu2", "mu3"), 1:unit.par))
      res.out$Sigma <- array(rbind(Sigma11.post, Sigma12.post, Sigma13.post,
                                   Sigma22.post, Sigma23.post, Sigma33.post),
                             dim=c(n.store, 6, unit.par),
                             dimnames=list(1:n.store, c("Sigma11",
                               "Sigma12", "Sigma13", "Sigma22", "Sigma23", "Sigma33"), 1:unit.par))
    }
    else {
      mu1.post <- matrix(res$pdSMu0, n.store, unit.par, byrow=TRUE)[,tmp$order.old]
      mu2.post <- matrix(res$pdSMu1, n.store, unit.par, byrow=TRUE)[,tmp$order.old]
      Sigma11.post <- matrix(res$pdSSig00, n.store, unit.par, byrow=TRUE)[,tmp$order.old]
      Sigma12.post <- matrix(res$pdSSig01, n.store, unit.par, byrow=TRUE)[,tmp$order.old]
      Sigma22.post <- matrix(res$pdSSig11, n.store, unit.par, byrow=TRUE)[,tmp$order.old]
      res.out$mu <- array(rbind(mu1.post, mu2.post), dim=c(n.store, 2, unit.par),
                          dimnames=list(1:n.store, c("mu1", "mu2"), 1:unit.par))
      res.out$Sigma <- array(rbind(Sigma11.post, Sigma12.post, Sigma22.post),
                             dim=c(n.store, 3, unit.par),
                             dimnames=list(1:n.store, c("Sigma11", "Sigma12", "Sigma22"), 1:unit.par))
    }
    if (alpha.update)
      res.out$alpha <- matrix(res$pdSa, n.store, unit.a, byrow=TRUE)
    else
      res.out$alpha <- alpha
    res.out$nstar <- matrix(res$pdSn, n.store, unit.a, byrow=TRUE)
  }

  if (context)
    class(res.out) <- c("ecoNPX", "ecoNP", "eco")
  else
      class(res.out) <- c("ecoNP", "eco")
  return(res.out)
}


