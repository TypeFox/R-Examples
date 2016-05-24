n4means <- function(delta, sigma, alpha=0.05, power=0.8, AR=1, two.tailed=TRUE, digits=3)
{
#Error Checking
if (sigma <=0) 
        stop("Sorry, the specified value of sigma must be strictly positive...")

if ((alpha >= 1) || (alpha <= 0) || (power <= 0) || (power >= 1))
        stop("Sorry, the alpha and power must lie within (0,1)")

if (AR <=0) 
        stop("Sorry, the specified value of the Allocation Ratio (AR) must be strictly positive...")

#Initialize Parameters
r <- NULL
r$delta <- delta; r$sigma <- sigma; r$alpha <- alpha; r$power <- power;
r$two.tailed <- two.tailed; r$AR <- AR;

#Calculate nE and nC for one or two tailed alpha
if (two.tailed)
{
r$n <- 2*sigma^2*(qnorm(1 - alpha/2) + qnorm(power))^2/(delta^2);
r$nE <- (1/2)*r$n*(1 + (1/AR));
r$nC <- (1/2)*r$n*(1 + AR);
}

if (!two.tailed)
{
r$n <- 2*sigma^2*(qnorm(1 - alpha) + qnorm(power))^2/(delta^2);
r$nE <- (1/2)*r$n*(1 + (1/AR));
r$nC <- (1/2)*r$n*(1 + AR);
}

class(r) <- "n4means";
return(r);

}

#Print method
print.n4means <- function(x, ...)
{
cat("The required sample size is a minimum of ", ceiling(x$nE), " individuals in the Experimental Group \n", sep="")
cat(" and a minimum of ", ceiling(x$nC), " individuals in the Control Group. \n", sep="")
}

#Summary Method
summary.n4means <- function(object, ...)
{
cat("Sample Size Calculation for Difference Between Means of Two Populations", "\n \n")
cat("Assuming:", "\n")
cat("Desired Minimum Detectable Difference between groups, delta = ", object$delta, "\n")
cat("Sigma = ", object$sigma, "\n");
cat("Type I Error Rate (alpha) = ", object$alpha, " and Power = ", object$power, "\n \n",sep="")

cat("The required sample size is a minimum of ", ceiling(object$nE), " individuals in the Experimental Group \n", sep="")
cat(" and a minimum of ", ceiling(object$nC), " individuals in the Control Group. \n", sep="")

}