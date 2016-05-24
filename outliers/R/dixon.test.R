"dixon.test" <-
function (x,type=0,opposite=FALSE,two.sided=TRUE)
{
DNAME <- deparse(substitute(x))
x <- sort(x[complete.cases(x)])

n <- length(x);

if ((type == 10 || type == 0) & (n < 3 || n > 30)) stop("Sample size must be in range 3-30");
if (type == 11 & (n < 4 || n > 30)) stop("Sample size must be in range 4-30");
if (type == 12 & (n < 5 || n > 30)) stop("Sample size must be in range 5-30");
if (type == 20 & (n < 4 || n > 30)) stop("Sample size must be in range 4-30");
if (type == 21 & (n < 5 || n > 30)) stop("Sample size must be in range 5-30");
if (type == 22 & (n < 6 || n > 30)) stop("Sample size must be in range 6-30");

if (sum(c(0,10,11,12,20,21,22) == type) == 0) stop ("Incorrect type");

if (type == 0) {
if (n <= 7) type <- 10
else if (n > 7 & n <=10) type <- 11
else if (n > 10 & n <=13) type <-21
else type <-22
} 

if (xor(((x[n] - mean(x)) < (mean(x) - x[1])),opposite)) {
alt = paste("lowest value",x[1],"is an outlier");
if ( type == 10 ) { Q = (x[2]-x[1])/(x[n]-x[1]) }
else
if ( type == 11 ) { Q = (x[2]-x[1])/(x[n-1]-x[1]) }
else
if ( type == 12 ) { Q = (x[2]-x[1])/(x[n-2]-x[1]) }
else
if ( type == 20 ) { Q = (x[3]-x[1])/(x[n]-x[1]) }
else
if ( type == 21 ) { Q = (x[3]-x[1])/(x[n-1]-x[1]) }
else
{ Q = (x[3]-x[1])/(x[n-2]-x[1]) }

}
else{
alt = paste("highest value",x[n],"is an outlier");
if ( type == 10 ) { Q = (x[n]-x[n-1])/(x[n]-x[1]) }
else
if ( type == 11 ) { Q = (x[n]-x[n-1])/(x[n]-x[2]) }
else
if ( type == 12 ) { Q = (x[n]-x[n-1])/(x[n]-x[3]) }
else
if ( type == 20 ) { Q = (x[n]-x[n-2])/(x[n]-x[1]) }
else
if ( type == 21 ) { Q = (x[n]-x[n-2])/(x[n]-x[2]) }
else
{ Q = (x[n]-x[n-2])/(x[n]-x[3]) }

}

pval <- pdixon(Q,n,type);

if (two.sided) {
pval <- 2* pval;
if (pval > 1) pval <- 2 - pval;
}

RVAL <- list(statistic = c(Q = Q), alternative = alt, p.value = pval, method = "Dixon test for outliers", 
    data.name = DNAME)
class(RVAL) <- "htest"
return(RVAL)
}

