##library(coda)

##takes args like MASS's mvrnorm
##rmvn <- function(n, mu=0, V = matrix(1))
##{
##  p <- length(mu)
##  if(any(is.na(match(dim(V),p))))
##    stop("Dimension problem!")
##  D <- chol(V)
##  matrix(rnorm(n*p), ncol=p) %*% D + rep(mu,rep(n,p))
##}

##takes mu and lower chol of cov matrix
##rmvn <- function(mu, chollt){
##  stdnorm <- as.matrix(rnorm(length(mu)))
##  chollt%*%stdnorm+as.matrix(mu)
##}


##final.sampler <- function (dataframe, totalsample, burnin, thinning)
##{
##select <- seq(from=burnin+1, to=totalsample, by=thinning)
##temp <- dataframe[select,]
##temp
##}

##gibbs.sample <- function(file, burnin, thinning)
##{
##dataframe <- read.table(file, header=F)
##select <- seq(from=burnin+1, to = nrow(dataframe), by=thinning)
##output <- dataframe[select,]
##output <- as.matrix(output)
##output
##}

##param.quantiles <- function(gibbsout)
##{
##
##  if (is.data.frame(gibbsout)) {
##    N <- ncol(gibbsout)
##    temp <- matrix(rep(0,times=N*3), ncol=3, byrow=T)
##
##	for(i in 1:N)
##	{
##	temp[i,] <- quantile(gibbsout[,i], c(0.50, 0.025, 0.975))
##      }
##  }
##  else {
##    N <- 1
##    temp <- quantile(gibbsout, c(0.50,0.025,0.975))
##  }
##  
##  temp
##}

##bayes.normal.reference <- function(sample.mean, sample.var, sample.size, NITER) {
##  n <- sample.size
##  y.bar <- sample.mean
##  s.sq <- sample.var
##  sigmasq.post <- 1/rgamma(NITER, (n-1)/2, (n-1)*s.sq/2)
##  mu.post <- rep(0,times=NITER)
##  for (i in 1:NITER) {
##    mu.post[i] <- rnorm(1, y.bar, sqrt(sigmasq.post[i]/n))
##  }
##  out <- cbind(mu.post, sigmasq.post)
##  out <- as.data.frame(out)
##  names(out) <- c("mu.posterior","sigmasq.posterior")
##  return(out)
##}


##bayesLMRef <- function(lm.obj, n.samples) {
##  lm.obj.data <- model.frame(lm.obj)
##  lm.obj.cov <- vcov(lm.obj)
##  D.inv = diag(diag(solve(lm.objcov)))
##  unit.corr = sqrt(D.inv)%*%unit.cov%*%sqrt(D.inv)
##  coeff.mean <- coefficients(lm.obj)

##  unit.coeff.post = rmvn(1000, mu=unit.coeff.mean, Sigma=unit.cov)
##  unit.coeff.post.qnt = param.quantiles(unit.coeff.post)

##  sigmasq <- (summary(lm.obj)$sigma)^2

##  N <- nrow(lm.obj.data)
##  p <- ncol(lm.obj.cov)
##  Simulate from marginal posterior [sigma^2 | y]
##  sigmasq.post <- 1/rgamma(n.samples, (N-p)/2, (N-p)*sigmasq/2)
  
##  Simulate from conditional posterior [beta | sigma^2, y]
##  coeff.post <- matrix(0, nrow=n.samples,ncol=p)
##  for (i in 1:n.samples) {
##    Sigma.coeff <- (sigmasq.post[i]/sigmasq)*lm.obj.cov
##    coeff.post[i,] <- rmvn(1, coeff.mean, Sigma.coeff)
##  }

##  Parameter.post <- as.data.frame(cbind(coeff.post, sigmasq.post))
##  colnames(Parameter.post) <- c(names(coefficients(lm.obj)),"sigma^2")
##  mcmc(Parameter.post)
##  Parameter.summary <- data.frame(param.quantiles(Parameter.post),row.names=c(names(coefficients(lm.obj)),"sigma^2")) 
##  names(Parameter.summary) <- c("50%","2.5%","97.5%")
##  Parameter.summary
##}
 





