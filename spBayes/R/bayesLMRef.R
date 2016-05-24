##takes args like MASS's mvrnorm
rmvn <- function(n, mu=0, V = matrix(1))
{
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  matrix(rnorm(n*p), ncol=p) %*% D + rep(mu,rep(n,p))
}


bayesLMRef <- function(lm.obj, n.samples, ...) {

  ####################################################
  ##Check for unused args
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }

  if(missing(lm.obj)){stop("error: lm.obj must be specified")}
  if(class(lm.obj) != "lm"){stop("error: lm.obj must be of class lm")}
  if(missing(n.samples)){stop("error: n.samples must be specified")} 
  
  lm.obj.data <- model.frame(lm.obj)
  lm.obj.cov <- vcov(lm.obj)
  ##D.inv = diag(diag(solve(lm.objcov)))
  ##unit.corr = sqrt(D.inv)%*%unit.cov%*%sqrt(D.inv)
  coeff.mean <- coefficients(lm.obj)

  ##unit.coeff.post = rmvn(1000, mu=unit.coeff.mean, Sigma=unit.cov)
  ##unit.coeff.post.qnt = param.quantiles(unit.coeff.post)

  sigmasq <- (summary(lm.obj)$sigma)^2

  N <- nrow(lm.obj.data)
  p <- ncol(lm.obj.cov)
  ## Simulate from marginal posterior [sigma^2 | y]
  sigmasq.post <- 1/rgamma(n.samples, (N-p)/2, (N-p)*sigmasq/2)
  
  ## Simulate from conditional posterior [beta | sigma^2, y]
  coeff.post <- matrix(0, nrow=n.samples,ncol=p)
  for (i in 1:n.samples) {
    Sigma.coeff <- (sigmasq.post[i]/sigmasq)*lm.obj.cov
    coeff.post[i,] <- rmvn(1, coeff.mean, Sigma.coeff)
  }

  Parameter.post <- cbind(coeff.post, sigmasq.post)
  colnames(Parameter.post) <- c(names(coefficients(lm.obj)),"tau.sq")
  ##mcmc(Parameter.post)
  ##Parameter.summary <- data.frame(param.quantiles(Parameter.post),row.names=c(names(coefficients(lm.obj)),"sigma.sq")) 
  ##names(Parameter.summary) <- c("50%","2.5%","97.5%")
  ##Parameter.summary


  ##
  ##return
  ##
  out <- list()
  out$p.beta.tauSq.samples <- mcmc(Parameter.post)
  out$X <- as.matrix(model.matrix(lm.obj))
  out$Y <- as.matrix(model.extract(model.frame(lm.obj), "response"))
  out$n.samples <- n.samples
  class(out) <- "bayesLMRef"
  out
}
 
