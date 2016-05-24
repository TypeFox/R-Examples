epiTTest <- function(X,Y,alpha=0.05, pooled=FALSE, digits=3)
{
#Error Checking
if (!is.vector(X))
	stop("Sorry, X and Y must be vectors...")

if (!is.vector(Y))
	stop("Sorry, X and Y must be vectors...")

if ((alpha >= 1) || (alpha <= 0))
        stop("Sorry, the alpha must lie within (0,1)")

#Initialize Object and Calculate simple statistics
r <- NULL;
r$digits <- digits;
r$alpha <- alpha;
r$pooled <- pooled;
r$nx <- length(X)
r$ny <- length(Y)
r$mean.x <- mean(X)
r$mean.y <- mean(Y)
r$s.x <- sd(X)
r$s.y <- sd(Y)
r$d <- r$mean.x - r$mean.y;

#Uses pooled estimates of variance
if (pooled)
{
r$df <- r$nx + r$ny -2;
r$s2p <- ((r$nx - 1)*var(X) +(r$ny -1)*var(Y))/r$df
r$CIL <- r$d - qt(1-alpha/2, r$df)*sqrt(r$s2p*((1/r$nx) + (1/r$ny)));
r$CIU <- r$d + qt(1-alpha/2, r$df)*sqrt(r$s2p*((1/r$nx) + (1/r$ny)));
r$TStat <- abs(r$d)/sqrt(r$s2p*((1/r$nx) + (1/r$ny)));
r$p.value <- 2* (1 - pt(r$TStat, r$df))
}

#Uses Satter. Approximation
if (!pooled)
{
r$df <- floor( (var(X)/r$nx + var(Y)/r$ny)^2/( (var(X)/r$nx)^2/(r$nx -1) + (var(Y)/r$ny)^2/(r$ny-1) ));
r$CIL <- r$d - qt(1-alpha/2, r$df)*sqrt(var(X)/r$nx + var(Y)/r$ny)
r$CIU <- r$d + qt(1-alpha/2, r$df)*sqrt(var(X)/r$nx + var(Y)/r$ny)
r$TStat <- abs(r$d)/sqrt(var(X)/r$nx + var(Y)/r$ny);
r$p.value <- 2*(1 - pt(r$TStat, r$df))
}

class(r) <- "epiTTest"
return(r);
}

#Print Method
print.epiTTest <- function(x, ...)
{
cat("Epidemiological T-Test Analysis", "\n \n")
cat("Number of Observations in Group I: ", x$nx, "\n", sep="")
cat("Number of Observations in Group II: ", x$ny, "\n \n", sep="")
cat("Sample Difference Between Means: ", round(x$d, digits=x$digits), "\n", sep="")
cat(100*(1-x$alpha), "% Confidence Limits for true mean difference: [", round(x$CIL,digits=x$digits), ", ", round(x$CIU,digits=x$digits),"]", "\n \n", sep="");

if (x$pooled)
{
cat("Note: The above Analysis uses a pooled estimate of the variance. \n")
}

if (!x$pooled)
{
cat("Note: The above Analysis uses the Satterthwaite estimate of the variance. \n")
}

}

#Summary Method
summary.epiTTest <- function(object, ...)
{
cat("Epidemiological T-Test Analysis", "\n \n")
cat("Number of Observations in Group I: ", object$nx, "\n", sep="")
cat("Sample Mean for Group I: ", round(object$mean.x, digits=object$digits), " with a sample standard deviation of ", round(object$s.x, digits=object$digits), "\n", sep="");
cat("Number of Observations in Group II: ", object$ny, "\n", sep="")
cat("Sample Mean for Group II: ", round(object$mean.y, digits=object$digits), " with a sample standard deviation of ", round(object$s.y, digits=object$digits), "\n", sep="");

if (object$pooled)
{
cat("Sample Difference Between Means: ", round(object$d, digits=object$digits), " with a pooled standard deviation of ", round(sqrt(object$s2p), digits=object$digits), "\n", sep="");
cat(100*(1-object$alpha), "% Confidence Limits for true mean difference: [", round(object$CIL,digits=object$digits), ", ", round(object$CIU,digits=object$digits),"]", "\n \n", sep="");
cat("T-Statistic for H0: difference = 0 vs. HA: difference != 0 is ", round(object$TStat,digits=object$digits), " with a p.value of ", round(object$p.value, digits=object$digits), "\n \n", sep="");
cat("Note: The above Analysis uses a pooled estimate of the variance. \n")

}

if (!object$pooled)
{
cat("Sample Difference Between Means: ", round(object$d, digits=object$digits), "\n", sep="")
cat(100*(1-object$alpha), "% Confidence Limits for true mean difference: [", round(object$CIL,digits=object$digits), ", ", round(object$CIU,digits=object$digits),"]", "\n \n", sep="");
cat("T-Statistic for H0: difference = 0 vs. HA: difference != 0 is ", round(object$TStat,digits=object$digits), " with a p.value of ", round(object$p.value, digits=object$digits), "\n \n", sep="");
cat("Note: The above Analysis uses the Satterthwaite estimate of the variance. \n")
}

}