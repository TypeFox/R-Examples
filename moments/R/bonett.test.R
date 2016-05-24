"bonett.test" <-
function (x, alternative=c("two.sided","less","greater"))
{
     DNAME <- deparse(substitute(x))
     x <- sort(x[complete.cases(x)])
     n <- length(x)
s <- match.arg(alternative)
alter <- switch(s, two.sided=0, less=1, greater=2)
rho <- sqrt(sum((x-mean(x))^2)/n);
tau <- sum(abs(x-mean(x)))/n;
omega <- 13.29*(log(rho)-log(tau));
z <- sqrt(n+2)*(omega-3)/3.54;
     pval <- pnorm(z, lower.tail = FALSE)
if (alter == 0) {
pval <- 2*pval
if (pval > 1) pval<-2-pval
alt <- "kurtosis is not equal to sqrt(2/pi)"
}
else if (alter == 1)
{
alt <- "kurtosis is greater than sqrt(2/pi)"
}
else
{
pval <- 1-pval
alt <- "kurtosis is lower than sqrt(2/pi)"
}
     RVAL <- list(statistic = c(tau = tau, z = z), alternative = alt, 
p.value = pval, method = "Bonett-Seier test for Geary kurtosis",
         data.name = DNAME)
     class(RVAL) <- "htest"
     return(RVAL)
}

