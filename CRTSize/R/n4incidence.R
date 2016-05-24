n4incidence <- function(le, lc, m, t, CV, alpha=0.05, power=0.8, AR=1, two.tailed=TRUE, digits=3)
{
#Error Checking
if ( (le <=0) || (lc <=0) )
        stop("Sorry, the specified value of le and lc must be strictly positive...")

if ((alpha >= 1) || (alpha <= 0) || (power <= 0) || (power >= 1))
        stop("Sorry, the alpha and power must lie within (0,1)...")

if (t <= 0)
        stop("Sorry, the value of t must lie be strictly positive...")

if (AR <=0) 
        stop("Sorry, the specified value of the Allocation Ratio (AR) must be strictly positive...")

if (m <=1) 
        stop("Sorry, the (average) cluster size, m, should be greater than one...")

#If m is a decimal, round up to generate a more conservative sample size.
m <- ceiling(m);

#Initialize Parameters
r <- NULL
r$le <- le; r$lc <- lc; r$t <- t; r$CV <- CV; r$m <- m; 
r$alpha <- alpha; r$power <- power; r$two.tailed <- two.tailed; r$AR <- AR;

#Calculate total number of clusters n,
if (two.tailed)
{
T <- t*m;
IFt <- 1 + ((CV^2)*(le^2 + lc^2)*T)/(le + lc);
r$n <- ( (qnorm(1 - alpha/2) + qnorm(power))^2*(le + lc))/(T*(le - lc)^2);
r$n <- IFt*r$n;
#Perform iterative sample size, using the T statistic for small n.

if (r$n < 30)
{
nTemp <- 0;
while (abs(r$n - nTemp) > 1)
{
nTemp <- r$n;
IFt <- 1 + ((CV^2)*(le^2 + lc^2)*T)/(le + lc);
r$n <- ( ( qt( (1 - alpha/2), df=( 2*(nTemp - 1))) + qt(power, df=(2*(nTemp - 1))) )^2*(le + lc))/(T*(le - lc)^2);
r$n <- IFt*r$n;
}

}

}

if (!two.tailed)
{
IFt <- 1 + ((CV^2)*(le^2 + lc^2)*T)/(le + lc);
r$n <- ( (qnorm(1 - alpha) + qnorm(power))^2*(le + lc))/(T*(le - lc)^2);
r$n <- IFt*r$n;

if (r$n < 30)
{

nTemp <- 0;
while (abs(r$n - nTemp) > 1)
{
nTemp <- r$n;
IFt <- 1 + ((CV^2)*(le^2 + lc^2)*T)/(le + lc);
r$n <- ( ( qt( (1 - alpha), df=( 2*(nTemp - 1))) + qt(power, df=(2*(nTemp - 1))) )^2*(le + lc))/(T*(le - lc)^2);
r$n <- IFt*r$n;
}

}

}

#Adjust for allocation ratio;
r$nE <- (1/2)*r$n*(1 + (1/AR));
r$nC <- (1/2)*r$n*(1 + AR);

class(r) <- "n4incidence";
return(r);

}

#Print method
print.n4incidence <- function(x, ...)
{
cat("The required sample size is a minimum of ", ceiling(x$nE), " clusters of size ", x$m, " in the Experimental Group \n", sep="")
cat(" and a minimum of ", ceiling(x$nC), " clusters (size ", x$m, ") in the Control Group, followed for time period of length ", x$t, "\n", sep="")
}

#Summary Method
summary.n4incidence <- function(object, ...)
{
cat("Sample Size Calculation for Comparison of Incidence Rates", "\n \n")
cat("Assuming:", "\n")
cat("Treatment Incidence Rate, le = ", object$le, "\n")
cat("Control Incidence Rate, lc = ", object$lc, "\n")
cat("Time Period, t = ", object$t, "\n")
cat("Cluster Size (average) = ", object$m, "\n");
cat("Coefficient of Variation, CV = ", object$CV, "\n")
cat("Type I Error Rate (alpha) = ", object$alpha, " and Power = ", object$power, "\n \n",sep="")

cat("The required sample size is a minimum of ", ceiling(object$nE), " clusters of size ", object$m, " in the Experimental Group \n", sep="")
cat(" and a minimum of ", ceiling(object$nC), " clusters (size ", object$m, ") in the Control Group.  Followed for time period of ", object$t, "\n", sep="")

}