######################################################################
# Function to compute MLE and backtransformed logit confidence
# interval for the pooled sample situation.
# In the simplest case a sensitivity and specificity of one
# is assumed, i.e. no misclassification occur.
#
# Warning: No adjusting is performed for the x=0 or the x=m
# situation.
#
# Params:
#   x - Number of positive pools (can be a vector)
#   k - pool size
#   n - number of pools
######################################################################

poolbinom.waldtype <- function(x, k, n, conf.level = 0.95,transformation=c("none","logit")) {

  #Handle possible vector arguments
  if ((length(x) != length(k)) | (length(x) != length(n))) {
    m <- cbind(x = x, k = k, n = n)
    x <- m[, "x"]
    k <- m[, "k"]
    n <- m[, "n"]
  }

  transformation <- match.arg(transformation,c("none","logit"))
  alpha <- 1-conf.level

  #Transformation function between observed Bin(m,theta) and pi of pooled
  g <- function(pi, x, k)  return(1 - (1- pi)^(1/k))

  #Pooled MLE
  pi.hat <- g(x/n,x=x, k=k)

  #Compute confidence interval and backtransform
  if (transformation == "logit") {
    #Look at where the logit transform can have problems
    xIs0 <- x == 0
    xIsn <- x == n
            
    #Delta-rule computed se of logit(pi.hat) -- computed by Maple
    logitpi.hat.se <- sqrt( x / k^2 / (n - x) / (-1 + ((n - x) / n)^(1 / k))^2 / n )

    logit.cis <- qlogis(pi.hat) + outer(qnorm(1-alpha/2) * logitpi.hat.se, c(-1,1))
    cis <- plogis(logit.cis)

    #Fix problematic intervals
    if (sum(xIs0)>0) {
      cis[xIs0,1] <- 0
#    cis[xIs0,2] <- 1 - g((alpha[xIs0]/2)^(1/n[xIs0]),x=x[xIs0],k=k[xIs0])
      cis[xIs0,2] <- poolbinom.lrt(x=x[xIs0], k=k[xIs0], n=n[xIs0])[,"upper"]
    }
    if (sum(xIsn)>0) {
#    cis[xIsn,1] <- g((alpha[xIsn]/2)^(1/n[xIsn]),x=x[xIsn],k=k[xIsn])
      cis[xIsn,1] <- poolbinom.lrt(x=x[xIsn], k=k[xIsn], n=n[xIsn])[,"lower"]
      cis[xIsn,2] <- 1
    }
  }
  #On original scale - no logit
  if (transformation == "none") {
    se <- sqrt((x/n)*(1-x/n)^(2/k-1)/k^2/n)
    cis <- pi.hat + outer(qnorm(1-alpha/2)*se, c(-1,1))
  }
  
  #Wrap up result in binom style
  res <- data.frame(method="logit",x=x,k=k, n=n, mle=pi.hat,lower=cis[,1],upper=cis[,2])
  return(res)
}

poolbinom.wald <- function(x, k, n, conf.level = 0.95) {
  poolbinom.waldtype(x=x, k=k, n=n, conf.level=conf.level,transformation=,"none")
}


poolbinom.logit <- function(x, k, n, conf.level = 0.95) {
  poolbinom.waldtype(x=x, k=k, n=n, conf.level=conf.level,transformation=,"logit")
}
