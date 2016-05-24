univar <- function(X,alpha=0.05, mu0 = 0, shapiro=FALSE, digits=3)
{
#Error Checking
if (!is.vector(X))
	stop("Sorry, X must be a vector...");

if ((alpha >= 1) || (alpha <= 0))
        stop("Sorry, the alpha must lie within (0,1)");
#Initialize parameters
Result <- NULL;
Result$mu0 <- mu0;
Result$n <- length(X);
Result$alpha <- alpha;
Result$digits <- digits;
Result$shapiro <- shapiro;

#Compute basic statistics
Result$mean <- mean(X);
Result$median <- median(X);
Result$min <- min(X);
Result$max <- max(X);
Result$s <- sd(X);
Result$var <- var(X)
Result$CIL <- Result$mean - qt(1-alpha/2,length(X)-1)*Result$s/sqrt(length(X));
Result$CIU <- Result$mean + qt(1-alpha/2,length(X)-1)*Result$s/sqrt(length(X));

#Hypothesis test for mu
Result$test <- sqrt(length(X))*abs(Result$mean - mu0)/Result$s;
Result$p.value <- 2*pt(Result$test, Result$n - 1, lower.tail=FALSE);

#Shapiro test if required
!(shapiro)
{
Result$shapiro.statistic <- shapiro.test(X)$statistic;
Result$shapiro.p.value <- shapiro.test(X)$p.value;
}
class(Result) <- "univar";
return(Result);
}

#Print Method
print.univar <- function(x, ...)
{
cat("\n","Univariate Summary", "\n \n", sep="")
cat("Sample Size: ", x$n, "\n", sep="" )
cat("Sample Mean:", round(x$mean,digits=x$digits), "\n", sep=" ") 
cat("Sample Median:", round(x$median,digits=x$digits), "\n", sep=" ") 
cat("Sample Standard Deviation:", round(x$s,digits=x$digits), "\n", sep=" ") 
}

#Summary Method
summary.univar <- function(object, ...)
{
cat("\n", "Basic Univariate Output", "\n \n", sep="")
cat("Sample Size: ", object$n, "\n", sep="" )
cat("Sample Mean:", round(object$mean,digits=object$digits), "\n", sep=" ") 
cat("Sample Median:", round(object$median,digits=object$digits), "\n", sep=" ") 
cat("The sample values range from:  ", round(object$min, digits=object$digits), " through ", round(object$max, digits=object$digits), "\n", sep="");
cat("Sample Standard Deviation:", round(object$s,digits=object$digits), "\n", sep=" ") 
cat("Sample Variance:", round(object$var,digits=object$digits), "\n \n ", sep=" ") 
cat(100*(1-object$alpha), "% Confidence Limits for mu are: [", round(object$CIL,digits=object$digits), ", ", round(object$CIU,digits=object$digits),"]", "\n", sep="");
cat("p-value for the Test of H0: mu = ", object$mu0, " vs. HA: mu != ", object$mu0, " is ", round(object$p.value,digits=object$digits), "\n \n", sep="");
cat("Note: P-Values and Confidence Intervals for mu assume an approximate normal distribution!", "\n", sep="")

if (object$shapiro)
cat("Shapiro-Wilks Statistic: ", round(object$shapiro.statistic, digits=object$digits), " with an approximate p-value of: ", round(object$shapiro.p.value,digits=object$digits), "\n", sep="");
}