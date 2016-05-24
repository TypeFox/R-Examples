jb.test <-
function(x)
{
n <- length(x)
avgx <- mean(x)
x.avgx <- x-avgx
x.avgx2 <- x.avgx^2
K <- n*sum(x.avgx^4)/(sum(x.avgx2)^2)
S <- (sum(x.avgx^3)/n)/(sum(x.avgx2)/n)^(3/2)
JB <- (n/6)*(S^2 + 0.25*((K-3)^2))
pval <- 1-pchisq(JB, df = 2)

out <- list()
out$skewness <- S
out$kurtosis <- K
out$statistic <- JB
out$df <- 2
out$p.value <- pval

return(out)
}
