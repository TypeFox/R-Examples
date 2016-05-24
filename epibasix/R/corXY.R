corXY <- function(X,Y, alpha=0.05, rho0 = 0, HA="not.equal", digits=3)
{
#Check for errors
if (length(X) != length(Y)) 
        stop("Sorry, this function Requires Two Vectors of the Same Length")

if ((alpha >= 1) || (alpha <= 0))
        stop("Sorry, the alpha must lie within (0,1)")

#Initialize object and parameters
result <- NULL;
result$digits <- digits;
result$alpha <- alpha;
result$rho0 <- rho0;
result$HA = HA;

result$rho <- cor(X,Y);
result$n <- length(X);

#Fisher's Z Transform and Hypothesis Testing Section
result$Z <- (1/2)*log((1+result$rho)/(1-result$rho));
result$Z0 <- (1/2)*log((1+rho0)/(1-rho0));

if (HA == "not.equal")
{
result$Test <- abs(result$Z - result$Z0)*sqrt(result$n-3);
result$p.value <- 2*(1 - pnorm(result$Test));
}

if (HA == "less.than")
{
result$Test <- abs(result$Z - result$Z0)*sqrt(result$n-3);
result$p.value <- pnorm(result$Test);
}

if (HA == "greater.than")
{
result$Test <- abs(result$Z - result$Z0)*sqrt(result$n-3);
result$p.value <- (1 - pnorm(result$Test));
}

#Compute Confidence Intervals
L <- result$Z - qnorm(1 - alpha/2)/sqrt(result$n-3);
U <- result$Z + qnorm(1 - alpha/2)/sqrt(result$n-3);

result$CIL <- (exp(2*L) -1)/(exp(2*L) + 1);
result$CIU <- (exp(2*U) -1)/(exp(2*U) + 1);

class(result) <- "corXY";

return(result);
}

#Print Method
print.corXY <- function(x, ...)
{
cat("\n","Correlation Summary", "\n \n", sep="")
cat("Sample Size: ", x$n, "\n", sep="" )
cat("Sample Correlation: ", round(x$rho, digits=x$digits), "\n", sep="")
cat(100*(1-x$alpha), "% Confidence Limits for rho are (using Fisher's Z Transformation) : [", round(x$CIL,digits=x$digits), ", ", round(x$CIU,digits=x$digits),"]", "\n", sep="");
}

#Summary Method
summary.corXY <- function(object, ...)
{
cat("\n","Correlation Details", "\n \n", sep="")
cat("Sample Size: ", object$n, "\n", sep="" )
cat("Sample Correlation: ", round(object$rho, digits=object$digits), "\n", sep="")
cat(100*(1-object$alpha), "% Confidence Limits for rho are (using Fisher's Z Transformation) : [", round(object$CIL,digits=object$digits), ", ", round(object$CIU,digits=object$digits),"]", "\n", sep="");
cat("Z statistic for H0: rho = ", object$rho0, " vs. HA: rho ", object$HA, " to ", object$rho0, ": ", round(object$Test, digits=object$digits), " , which has a p.value of ", round(object$p.value, digits=object$digits), "\n", sep="");
}