fixedMetaAnalRROR <- function(data, alpha=0.05)
{
if (!is.matrix(data))
	stop("Sorry data must be a matrix of OR/RR, 95 % Lower and Upper Limits from Previous Studies")

if (ncol(data) != 3)
	stop("Data must have 3 columns, Odds Ratio/Relative Risk, Lower Limit and Upper Limit from Previous Studies")

if ((alpha >= 1) || (alpha <= 0))
        stop("Sorry, the alpha must lie within (0,1)")

X <- NULL;
X$data <- data; X$alpha <- alpha;
logRR <- log(data[,1])
logL <- log(data[,2])
logU <- log(data[,3])


colnames(X$data) <- c("OR/Rel Risk", "Lower Limit", "Upper Limit");

selogRR <- (logU - logRR)/1.96;
varlogRR <- selogRR^2;
w <- 1/varlogRR;

Z <- -qnorm(alpha/2)
X$thetaF <- sum(logRR*w)/sum(w);
X$uF <- X$thetaF + Z/sqrt(sum(w));
X$lF <- X$thetaF - Z/sqrt(sum(w));
X$Var <- 1/sum(w);

if ( (X$uF < 0) && (X$lF < 0) || (X$uF > 0) && (X$lF > 0) )
{
X$Sig <- 1;
}
else
{
X$Sig <- 0;
}


class(X) <- "fixedMetaAnalRROR";
return(X);
}


#Print method
print.fixedMetaAnalRROR <- function(x, ...)
{
cat("The fixed effects Relative Risk/Odds Ratio is ", exp(x$thetaF), "\n", sep="");
cat("With ", (1 - x$alpha)*100,  "% Confidence Limits: (", exp(x$lF), ", ", exp(x$uF),")", sep="");
}

#Summary Method
summary.fixedMetaAnalRROR <- function(object, ...)
{
cat("Original Data Matrix:", "\n \n", sep="")
print(object$data)
cat(" \n The fixed effects Relative Risk/Odds Ratio is ", exp(object$thetaF), "\n", sep="");
cat("With ", (1 - object$alpha)*100,  "% Confidence Limits: (", exp(object$lF), ", ", exp(object$uF), ") \n \n",sep="");
if (object$Sig == 1)
{
cat("This result is statistically significant at the 5 percent level.", sep="");
}
if (object$Sig == 0)
{
cat("This result is not statistically significant at the 5 percent level.", sep="");
}
}