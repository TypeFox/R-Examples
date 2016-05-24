###
### Bayesian probit with nonignorable missing outcomes and
### multi-valued treatments
###

NIbprobit <- function(formula, Xo, Xr, data = parent.frame(),
                      n.draws = 5000, insample = FALSE,
                      param = TRUE, mda = TRUE,
                      p.mean.o = 0, p.prec.o = 0.01,
                      p.mean.r = 0, p.prec.r = 0.01,
                      coef.start.o = 0, coef.start.r = 0,
                      burnin = 0, thin = 0, verbose = TRUE) {  

  ## getting Y and D
  call <- match.call()
  tm <- terms(formula)
  attr(tm, "intercept") <- 0
  mf <- model.frame(tm, data = data, na.action = 'na.pass')
  D <- model.matrix(tm, data = mf)
  if (max(D) > 1 || min(D) < 0)
    stop("the treatment variable should be a factor variable.")
  Y <- model.response(mf)
  m <- ncol(D) # number of treatment levels including control
  ## getting Xo and Xr
  tm <- terms(Xo)
  attr(tm, "intercept") <- 1
  Xo <- model.matrix(tm, data = data, na.action = 'na.pass')
  Xo <- Xo[,(colnames(Xo) != "(Intercept)")]
  tm <- terms(Xr)
  attr(tm, "intercept") <- 1
  Xr <- model.matrix(tm, data = data, na.action = 'na.pass')
  Xr <- Xr[,(colnames(Xr) != "(Intercept)")]
  ## taking care of NA's in D and X
  ind <- complete.cases(cbind(D, Xo, Xr))
  Y <- Y[ind]
  D <- D[ind,]
  Xo <- Xo[ind,]
  Xr <- Xr[ind,]
  R <- (!is.na(Y))*1
  Y[is.na(Y)] <- rbinom(sum(is.na(Y)), size = 1, prob = 0.5)
  cnameso <- c(colnames(D), colnames(Xo))
  cnamesr <- c("1-Y", "Y", colnames(Xr))
  Xo <- cbind(D, Xo)
  colnames(Xo) <- cnameso
  Xr <- cbind(1-Y, Y, Xr)
  colnames(Xr) <- cnamesr
  
  res <- list(call = call, Y = Y, Xo = Xo, Xr = Xr, n.draws = n.draws)

  n <- length(Y)
  ncovo <- ncol(Xo)
  ncovr <- ncol(Xr)
  ## starting values
  if(length(coef.start.o) != ncovo)
    coef.start.o <- rep(coef.start.o, ncovo)
  if(length(coef.start.r) != ncovr)
    coef.start.r <- rep(coef.start.r, ncovr)
 
  ## prior
  if(length(p.mean.o) != ncovo)
    p.mean.o <- rep(p.mean.o, ncovo)
  if(length(p.mean.r) != ncovr)
    p.mean.r <- rep(p.mean.r, ncovr)
  if(!is.matrix(p.prec.o))
    p.prec.o <- diag(p.prec.o, ncovo)
  if(!is.matrix(p.prec.r))
    p.prec.r <- diag(p.prec.r, ncovr)
  
  ## checking thinnig and burnin intervals
  if (n.draws <= 0)
    stop("`n.draws' should be a positive integer.")
  if (burnin < 0 || burnin >= n.draws)
    stop("`burnin' should be a non-negative integer less than `n.draws'.")
  if (thin < 0 || thin >= n.draws)
    stop("`thin' should be a non-negative integer less than `n.draws'.")
  keep <- thin + 1

  ## calling C function to do MCMC
  par <- .C("NIbprobit",
            as.integer(Y), as.integer(R), 
            as.double(Xo), as.double(Xr),
            as.double(coef.start.o), as.double(coef.start.r),
            as.integer(n), as.integer(ncovo), as.integer(ncovr),
            as.integer(m), as.double(p.mean.o), as.double(p.mean.r),
            as.double(p.prec.o), as.double(p.prec.r), 
            as.integer(insample), as.integer(param), as.integer(mda),
            as.integer(n.draws), as.integer(burnin),
            as.integer(keep), as.integer(verbose),
            coef.o = double(ncovo*(ceiling((n.draws-burnin)/keep))),
            coef.r = double(ncovr*(ceiling((n.draws-burnin)/keep))),
            ATE = double((m-1)*(ceiling((n.draws-burnin)/keep))),
            BASE = double(m*(ceiling((n.draws-burnin)/keep))),
            PACKAGE="experiment")
  if (param) {
    res$coef.o <- matrix(par$coef.o, byrow = TRUE, ncol = ncovo)
    colnames(res$coef.o) <- colnames(Xo)
    res$coef.r <- matrix(par$coef.r, byrow = TRUE, ncol = ncovr)
    colnames(res$coef.r) <- colnames(Xr)
  }
  res$ATE <- matrix(par$ATE, byrow = TRUE, ncol = m-1)
  res$base <- matrix(par$BASE, byrow = TRUE, ncol = m)
  colnames(res$base) <- colnames(D)
  
  class(res) <- "NIbprobit"
  return(res)
}
