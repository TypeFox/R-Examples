"chisq.out.test" <-
function (x,variance=var(x),opposite=FALSE) 
{
DNAME <- deparse(substitute(x))
x <- sort(x[complete.cases(x)])
n <- length(x);
if (xor(((x[n] - mean(x)) < (mean(x) - x[1])),opposite)) {
alt = paste("lowest value",x[1],"is an outlier");
chi <- (mean(x)-x[1])^2/variance
}
else{
alt = paste("highest value",x[n],"is an outlier");
chi <- (x[n]-mean(x))^2/variance
}
pval <- 1-pchisq(chi,1)
RVAL <- list(statistic = c("X-squared" = chi), alternative = alt, p.value = pval, method = "chi-squared test for outlier", 
    data.name = DNAME)
class(RVAL) <- "htest"
return(RVAL)
}

