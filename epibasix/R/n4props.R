n4props <- function(pe,pc, alpha=0.05, power=0.8, AR=1, two.tailed=TRUE, digits=3)
{
#Error Checking
if ((pe >= 1) || (pe <= 0) || (pc <= 0) || (pe >= 1))
        stop("Sorry, the prior proportions must lie within (0,1)")

if ((alpha >= 1) || (alpha <= 0) || (power <= 0) || (power >= 1))
        stop("Sorry, the alpha and power must lie within (0,1)")

if (AR <=0) 
        stop("Sorry, the specified value of the Allocation Ratio (AR) must be strictly positive...")

#Initialize Parameters
r <- NULL;
r$pe <- pe; r$pc <- pc; r$digits <- digits;
r$alpha <- alpha; r$power <- power;
r$AR <- AR; r$two.tailed <- two.tailed; 

#One or two-tailed tests
if (two.tailed)
{
r$n <- (qnorm(1 - alpha/2) + qnorm(power))^2*(pe*(1-pe) + pc*(1-pc))/(pe - pc)^2;
r$nE = (1/2)*r$n*(1 + (1/AR));
r$nC = (1/2)*r$n*(1 + AR);
}

else
{
r$n <- (qnorm(1 - alpha) + qnorm(power))^2*(pe*(1-pe) + pc*(1-pc))/(pe - pc)^2;
r$nE = (1/2)*r$n*(1 + (1/AR));
r$nC = (1/2)*r$n*(1 + AR);
}

class(r) <- "n4props";
return(r);

}

#Print Method
print.n4props <- function(x, ...)
{
cat("The required sample size is a minimum of ", ceiling(x$nE), " individuals in the Experimental Group \n", sep="")
cat(" and a minimum of ", ceiling(x$nC), " individuals in the Control Group. \n", sep="")
}

#Summary Method
summary.n4props <- function(object, ...)
{
cat("Sample Size Calculation for Binary Outcomes", "\n \n")
cat("Assuming:", "\n")
cat("Proportion with Outcome in Experimental Group: ", object$pe, "\n")
cat("Proportion with Outcome in Control Group: ", object$pc, "\n")
cat("Type I Error Rate (alpha) = ", object$alpha, " and Power = ", object$power, "\n \n",sep="")

cat("The required sample size is a minimum of ", ceiling(object$nE), " individuals in the Experimental Group \n", sep="")
cat(" and a minimum of ", ceiling(object$nC), " individuals in the Control Group. \n", sep="")

}