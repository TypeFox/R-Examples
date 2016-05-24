"grubbs.test" <-
function (x,type=10,opposite=FALSE,two.sided=FALSE)
{

if (sum(c(10,11,20) == type) == 0) stop ("Incorrect type");


DNAME <- deparse(substitute(x))
x <- sort(x[complete.cases(x)])

n <- length(x);

if (type == 11) 
{
g <- (x[n] - x[1])/sd(x);
u <- var(x[2:(n-1)])/var(x)*(n-3)/(n-1)

pval = 1-pgrubbs(g,n,type=11);

method <- "Grubbs test for two opposite outliers"

alt = paste(x[1],"and",x[n],"are outliers")

}

else if (type == 10)

{
if (xor(((x[n] - mean(x)) < (mean(x) - x[1])),opposite)) {
alt = paste("lowest value",x[1],"is an outlier");
o <- x[1];
d <- x[2:n];
}
else{
alt = paste("highest value",x[n],"is an outlier");
o <- x[n];
d <- x[1:(n-1)];

}

g <- abs(o - mean(x))/sd(x)
u <- var(d)/var(x)*(n-2)/(n-1)

pval <- 1-pgrubbs(g,n,type=10);

method <- "Grubbs test for one outlier"

}

else

{

if (xor(((x[n] - mean(x)) < (mean(x) - x[1])),opposite)) {
alt = paste("lowest values",x[1],",",x[2],"are outliers");
u <- var(x[3:n])/var(x)*(n-3)/(n-1);
}
else{
alt = paste("highest values",x[n-1],",",x[n],"are outliers");
u <- var(x[1:(n-2)])/var(x)*(n-3)/(n-1)
}

g <- NULL

pval <- pgrubbs(u,n,type=20);

method <- "Grubbs test for two outliers"

}

if (two.sided) {
        pval <- 2* pval;
        if (pval > 1) pval <- 2 - pval;
        }


RVAL <- list(statistic = c(G = g, U= u),
alternative = alt, p.value = pval, method = method, 
    data.name = DNAME)
class(RVAL) <- "htest"
return(RVAL)


}

