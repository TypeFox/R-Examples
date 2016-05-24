geno <-c(83,447,470)
loglik <- function(theta,x) {
    probs <- c(theta^2, 2*theta*(1-theta), (1-theta)^2)
    if (any (probs <=0)) return (Inf)
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
}
