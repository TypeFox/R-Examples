eco <- function(formula, data = parent.frame(), N = NULL, supplement = NULL,
                context = FALSE, mu0 = 0, tau0 = 2, nu0 = 4, S0 = 10,
                mu.start = 0, Sigma.start = 10, parameter = TRUE,
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
  if (length(mu.start)==1)
    mu.start <- rep(mu.start, ndim)
  else if (length(mu.start)!=ndim)
    stop("invalid inputs for mu.start")
  if (is.matrix(Sigma.start)) {
    if (any(dim(Sigma.start)!=ndim))
      stop("invalid inputs for Sigma.start")
  }
  else
    Sigma.start <- diag(Sigma.start, ndim)
  
  ## getting X, Y, and N
  mf <- match.call()
  tt <- terms(formula)
  attr(tt, "intercept") <- 0
  if (is.matrix(eval.parent(mf$data)))
    data <- as.data.frame(data)
  X <- model.matrix(tt, data)
  Y <- model.response(model.frame(tt, data = data))
  N <- eval(mf$N, data)
  
  # check data and modify inputs 
  tmp <- checkdata(X,Y, supplement, ndim)  
  bdd <- ecoBD(formula=formula, data=data)
  W1min <- bdd$Wmin[order(tmp$order.old)[1:nrow(tmp$d)],1,1]
  W1max <- bdd$Wmax[order(tmp$order.old)[1:nrow(tmp$d)],1,1]
 

  ## fitting the model
  n.store <- floor((n.draws-burnin)/(thin+1))
  unit.par <- 1
  unit.w <- tmp$n.samp+tmp$samp.X1+tmp$samp.X0 	
  n.w <- n.store * unit.w

  if (context) 
    res <- .C("cBaseecoX", as.double(tmp$d), as.integer(tmp$n.samp),
              as.integer(n.draws), as.integer(burnin), as.integer(thin+1),
              as.integer(verbose), as.integer(nu0), as.double(tau0),
              as.double(mu0), as.double(S0), as.double(mu.start),
              as.double(Sigma.start), as.integer(tmp$survey.yes),
              as.integer(tmp$survey.samp), as.double(tmp$survey.data),
              as.integer(tmp$X1type), as.integer(tmp$samp.X1),
              as.double(tmp$X1.W1), as.integer(tmp$X0type),
              as.integer(tmp$samp.X0), as.double(tmp$X0.W2),
              as.double(W1min), as.double(W1max),
              as.integer(parameter), as.integer(grid), 
              pdSMu0 = double(n.store), pdSMu1 = double(n.store), pdSMu2 = double(n.store),
              pdSSig00=double(n.store), pdSSig01=double(n.store), pdSSig02=double(n.store),
              pdSSig11=double(n.store), pdSSig12=double(n.store), pdSSig22=double(n.store),
              pdSW1=double(n.w), pdSW2=double(n.w), PACKAGE="eco")
  else 
    res <- .C("cBaseeco", as.double(tmp$d), as.integer(tmp$n.samp),
              as.integer(n.draws), as.integer(burnin), as.integer(thin+1),
              as.integer(verbose), as.integer(nu0), as.double(tau0),
              as.double(mu0), as.double(S0), as.double(mu.start),
              as.double(Sigma.start), as.integer(tmp$survey.yes),
              as.integer(tmp$survey.samp), as.double(tmp$survey.data),
              as.integer(tmp$X1type), as.integer(tmp$samp.X1),
              as.double(tmp$X1.W1), as.integer(tmp$X0type),
              as.integer(tmp$samp.X0), as.double(tmp$X0.W2),
              as.double(W1min), as.double(W1max),
              as.integer(parameter), as.integer(grid), 
              pdSMu0=double(n.store), pdSMu1=double(n.store), 
	      pdSSig00=double(n.store),
              pdSSig01=double(n.store), pdSSig11=double(n.store),
              pdSW1=double(n.w), pdSW2=double(n.w),
              PACKAGE="eco")
    
  W1.post <- matrix(res$pdSW1, n.store, unit.w, byrow=TRUE)[,tmp$order.old]
  W2.post <- matrix(res$pdSW2, n.store, unit.w, byrow=TRUE)[,tmp$order.old]
  W <- array(rbind(W1.post, W2.post), c(n.store, 2, unit.w))
  colnames(W) <- c("W1", "W2")
  res.out <- list(call = mf, X = X, Y = Y, N = N, W = W,
                  Wmin=bdd$Wmin[,1,], Wmax = bdd$Wmax[,1,],
                  burin = burnin, thin = thin, nu0 = nu0,
                  tau0 = tau0, mu0 = mu0, S0 = S0)
  
  if (parameter) 
    if (context) {
      res.out$mu <- cbind(matrix(res$pdSMu0, n.store, unit.par, byrow=TRUE),
                          matrix(res$pdSMu1, n.store, unit.par, byrow=TRUE), 
                          matrix(res$pdSMu2, n.store, unit.par, byrow=TRUE)) 
      colnames(res.out$mu) <- c("mu1", "mu2", "mu3")
      res.out$Sigma <- cbind(matrix(res$pdSSig00, n.store, unit.par, byrow=TRUE), 
                             matrix(res$pdSSig01, n.store, unit.par, byrow=TRUE),
                             matrix(res$pdSSig02, n.store, unit.par, byrow=TRUE),
                             matrix(res$pdSSig11, n.store, unit.par, byrow=TRUE),
                             matrix(res$pdSSig12, n.store, unit.par, byrow=TRUE),
                             matrix(res$pdSSig22, n.store, unit.par, byrow=TRUE))
      colnames(res.out$Sigma) <- c("Sigma11", "Sigma12", "Sigma13",
                                   "Sigma22", "Sigma23", "Sigma33")      
    }
    else {
      res.out$mu <- cbind(matrix(res$pdSMu0, n.store, unit.par, byrow=TRUE),
                          matrix(res$pdSMu1, n.store, unit.par, byrow=TRUE)) 
      colnames(res.out$mu) <- c("mu1", "mu2")
      res.out$Sigma <- cbind(matrix(res$pdSSig00, n.store, unit.par, byrow=TRUE), 
                             matrix(res$pdSSig01, n.store, unit.par, byrow=TRUE),
                             matrix(res$pdSSig11, n.store, unit.par, byrow=TRUE))
      colnames(res.out$Sigma) <- c("Sigma11", "Sigma12", "Sigma22")
    }

  if (context)
    class(res.out) <- c("ecoX","eco")
  else
    class(res.out) <- c("eco")
  
  return(res.out)

}


