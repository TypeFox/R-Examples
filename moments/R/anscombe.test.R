"anscombe.test" <-
function (x, alternative=c("two.sided","less","greater"))
{
     DNAME <- deparse(substitute(x))
     x <- sort(x[complete.cases(x)])
     n <- length(x)
s <- match.arg(alternative)
alter <- switch(s, two.sided=0, less=1, greater=2)
b <- n*sum( (x-mean(x))^4 )/(sum( (x-mean(x))^2 )^2);
eb2 <- 3*(n-1)/(n+1);
vb2 <- 24*n*(n-2)*(n-3)/ ((n+1)^2*(n+3)*(n+5));
m3 <- (6*(n^2-5*n+2)/((n+7)*(n+9)))*sqrt((6*(n+3)*(n+5))/(n*(n-2)*(n-3)));
a <- 6 + (8/m3) * (2/m3 + sqrt(1 + 4/m3^2));
xx <- (b-eb2)/sqrt(vb2);
z <- ( 1-2/(9*a)-( (1-2/a) / (1+xx*sqrt(2/(a-4))) )^(1/3))/ sqrt(2/(9*a));
     pval <- pnorm(z, lower.tail = FALSE)
if (alter == 0) {
pval <- 2*pval
if (pval > 1) pval<-2-pval
alt <- "kurtosis is not equal to 3"
}
else if (alter == 1)
{
alt <- "kurtosis is greater than 3"
}
else
{
pval <- 1-pval
alt <- "kurtosis is lower than 3"
}
     RVAL <- list(statistic = c(kurt = b, z = z), p.value = pval, 
alternative = alt, method = "Anscombe-Glynn kurtosis test",
         data.name = DNAME)
     class(RVAL) <- "htest"
     return(RVAL)
}

