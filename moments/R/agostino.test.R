"agostino.test" <-
function (x, alternative=c("two.sided","less","greater"))
{
     DNAME <- deparse(substitute(x))
     x <- sort(x[complete.cases(x)])
     n <- length(x)
s <- match.arg(alternative)
alter <- switch(s, two.sided=0, less=1, greater=2)
     if ((n < 8 || n > 46340))
         stop("sample size must be between 8 and 46340")
s3 <- (sum((x-mean(x))^3)/n)/(sum((x-mean(x))^2)/n)^(3/2)
y <- s3*sqrt((n+1)*(n+3)/(6*(n-2)))
b2 <- 3*(n*n+27*n-70)*(n+1)*(n+3)/((n-2)*(n+5)*(n+7)*(n+9))
w <- sqrt(-1+sqrt(2*(b2-1)));
d <- 1/sqrt(log(w));
a <- sqrt(2/(w*w-1));
z <- d*log(y/a+sqrt((y/a)^2+1));
     pval <- pnorm(z, lower.tail = FALSE)
if (alter == 0) {
pval <- 2*pval
if (pval > 1) pval<-2-pval
alt <- "data have a skewness"
}
else if (alter == 1)
{
alt <- "data have positive skewness"
}
else
{
pval <- 1-pval
alt <- "data have negative skewness"
}
     RVAL <- list(statistic = c(skew = s3, z = z), p.value = pval, 
alternative = alt, method = "D'Agostino skewness test",
         data.name = DNAME)
     class(RVAL) <- "htest"
     return(RVAL)
}

