counts <- c(1997,906,904,32)
# computes model probabilities from value of theta
theta2probs <- function(theta) {
    c(0.25 * (2 + theta), 
      0.25 * (1-theta), 
      0.25 * (1-theta), 
      0.25 * theta )
}

# computes log( L(theta | x) )
loglik1 <- function(theta,x) {
     (   x[1]       * log(0.25 * (2 + theta)) 
      + (x[2]+x[3]) * log(0.25 * (1-theta))
      +  x[4]       * log(0.25 * theta)
      )
}

mle <- nlmax(loglik1,p=0.5,x=counts); summary(mle)
theta.hat <- mle$estimate

loglik <- function(theta,x) {
    if (theta < 0 || theta > 1) { return (Inf) }
    dmultinom( x, size = sum(x),prob=theta2probs(theta), log=T )
}

mle <- nlmax(loglik,p=0.5,x=counts); summary(mle)
theta.hat <- mle$estimate

# test a specific value of theta vs best possible theta
testTheta <- function(theta,x) {
    chisqStat <- 2 * (loglik(theta.hat,x) - loglik(theta,x))
    chisqStat1 <- 2 * (loglik1(theta.hat,x) - loglik1(theta,x))
    p.value <- 1 - pchisq(chisqStat,df=1)
    return( c(statistic = chisqStat, statistic1 = chisqStat1, 
              p.value = p.value))
    }
chisq.test(counts,p=theta2probs(mle$estimate))
# so we can grab that statistic and redo the p-value:
X <- chisq.test(counts,p=theta2probs(mle$estimate))$statistic
1 - pchisq(X,df=3-1)  # df=3 for multinomial, 1 for model based on theta

# alternatively, we can do this manually:
o <- counts
e <- theta2probs(theta.hat) * sum(o)
testStats <- c(lrt = 2 * sum( o * log (o/e)), pearson= sum( (o-e)^2/e) )
testStats
1-pchisq(testStats,df=3-1)
