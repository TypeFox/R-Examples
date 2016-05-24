n4means <- function(delta, sigma, m, ICC, alpha=0.05, power=0.8, AR=1, two.tailed=TRUE, digits=3)
{
#Error Checking
if (sigma <=0) 
        stop("Sorry, the specified value of sigma must be strictly positive...")

if ((alpha >= 1) || (alpha <= 0) || (power <= 0) || (power >= 1))
        stop("Sorry, the alpha and power must lie within (0,1)")

if (ICC <= 0)
        stop("Sorry, the ICC must lie within (0,1)")

if (AR <=0) 
        stop("Sorry, the specified value of the Allocation Ratio (AR) must be strictly positive...")

if (m <=1) 
        stop("Sorry, the (average) cluster size, m, should be greater than one...")

#If m is a decimal, round up to generate a more conservative sample size.
m <- ceiling(m);

#Initialize Parameters
r <- NULL
r$delta <- delta; r$sigma <- sigma; r$ICC <- ICC; r$m <- m; 
r$alpha <- alpha; r$power <- power; r$two.tailed <- two.tailed; r$AR <- AR;

#Calculate total number of clusters n,
if (two.tailed)
{
r$n <- 2*sigma^2*(1 + (m - 1)*ICC)*(qnorm(1 - alpha/2) + qnorm(power))^2/(m*(delta^2));

#Perform iterative sample size, using the T statistic for small n.

if (r$n < 30)
{
nTemp <- 0;
while (abs(r$n - nTemp) > 1)
{
nTemp <- r$n;
r$n <- 2*sigma^2*(1 + (m - 1)*ICC)*((qt((1 - alpha/2), df=( 2*(nTemp - 1) ))) + qt(power, df=(2*(nTemp - 1))))^2/(m*(delta^2));
}

}

}

if (!two.tailed)
{
r$n <- 2*sigma^2*(1 + (m - 1)*ICC)*(qnorm(1 - alpha) + qnorm(power))^2/(m*(delta^2));

if (r$n < 30)
{

nTemp <- 0;
while (abs(r$n - nTemp) > 1)
{
nTemp <- r$n;
r$n <- 2*sigma^2*(1 + (m - 1)*ICC)*((qt((1 - alpha/2), df=( 2*(nTemp - 1) ))) + qt(power, df=(2*(nTemp - 1))))^2/(m*(delta^2));
}

}

}

#Adjust for allocation ratio;
r$nE <- (1/2)*r$n*(1 + (1/AR));
r$nC <- (1/2)*r$n*(1 + AR);

class(r) <- "n4means";
return(r);

}

#Print method
print.n4means <- function(x, ...)
{
cat("The required sample size is a minimum of ", ceiling(x$nE), " clusters of size ", x$m, " in the Experimental Group \n", sep="")
cat(" and a minimum of ", ceiling(x$nC), " clusters (size ", x$m, ") in the Control Group. \n", sep="")
}

#Summary Method
summary.n4means <- function(object, ...)
{
cat("Sample Size Calculation for Difference Between Means of Two Populations", "\n \n")
cat("Assuming:", "\n")
cat("Desired Minimum Detectable Difference between groups, delta = ", object$delta, "\n")
cat("Sigma = ", object$sigma, "\n");
cat("Cluster Size (average) = ", object$m, "\n");
cat("ICC = ", object$ICC, "\n");
cat("Type I Error Rate (alpha) = ", object$alpha, " and Power = ", object$power, "\n \n",sep="")

cat("The required sample size is a minimum of ", ceiling(object$nE), " clusters of size ", object$m, " in the Experimental Group \n", sep="")
cat(" and a minimum of ", ceiling(object$nC), " clusters (size ", object$m, ") in the Control Group. \n", sep="")

}