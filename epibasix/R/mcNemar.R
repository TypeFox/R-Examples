mcNemar <- function(X, alpha = 0.05, force=FALSE, digits=3)
{
#Error Checking
if (nrow(X) != ncol(X)) 
        stop("Sorry, 'X' must be a 2x2 Matrix...")

#For notation...
a <- X[1,1];
b <- X[1,2];
c <- X[2,1];
d <- X[2,2]
n <- sum(X);

if ((alpha >= 1) || (alpha <= 0))
        stop("Sorry, the alpha must lie within (0,1)")

#To proceed with function with insufficient data, set force=TRUE
if (!force)
{
if (b+c <= 30)
	stop("The number of discordant pairs, must be greater than 30 to ensure validity.")   
}

#Names for Rows and Columns
if ( (is.null(rownames(X))) && (is.null(colnames(X))) )
{
colnames(X) <- c("Exposed Person: Disease Present", "Exposed Person: Disease Absent")
rownames(X) <- c("Control Person: Disease Present", "Control Person: Disease Absent");
}

#Initialize Parameters
result <- NULL;

result$X <- X;
result$alpha <- alpha;
result$digits <- digits;

#McNemar's Odds Ratio
result$ORMc <- b/c;

p <- b/(b+c);
pL <- p - 1/(2*(b+c)) - qnorm(1-alpha/2)*sqrt((b*c)/((b+c)^3)); 
pU <- p + 1/(2*(b+c)) + qnorm(1-alpha/2)*sqrt((b*c)/((b+c)^3));

#OR Confidence Limits
result$ORMc.CIL <- pL/(1-pL);
result$ORMc.CIU <- pU/(1-pU);

#Risk Difference and CI
result$rd <- (b - c)/n;
result$rd.CIL <- result$rd - (1/n) - qnorm(1-alpha/2)*(sqrt((a + d)*(b+c)+4*b*c)/(n*sqrt(n)));
result$rd.CIU <- result$rd + (1/n) + qnorm(1-alpha/2)*(sqrt((a + d)*(b+c)+4*b*c)/(n*sqrt(n)));

#McNemarStat, with Continuity Correction
result$XMc <- ((abs(b-c) -1)^2)/(b+c);
result$XMc.p.value <- 1 - pchisq(result$XMc,1);   

class(result) <- "mcNemar";
return(result)
}

#Print Method
print.mcNemar <- function(x, ...)
{
cat("\n","Matched Pairs Analysis: McNemar's Chi^2 Statistic and Odds Ratio", "\n \n", sep="")
cat("McNemar's Chi^2 Statistic (corrected for continuity) = ", round(x$XMc, digits=x$digits), " which has a p-value of: ", round(x$XMc.p.value, digits=x$digits), "\n \n",sep="")
cat("McNemar's Odds Ratio (b/c): ", round(x$ORMc, digits=x$digits), "\n", sep="") 
cat(100*(1-x$alpha), "% Confidence Limits for the OR are: [", round(x$ORMc.CIL,digits=x$digits), ", ", round(x$ORMc.CIU,digits=x$digits),"]", "\n", sep="");

}

#Summary Method
summary.mcNemar <- function(object, ...)
{

cat("\n", "Matched Pairs Analysis: McNemar's Statistic and Odds Ratio (Detailed Summary):", "\n \n", sep="")
print(object$X);
cat("\n");
cat("Entries in above matrix correspond to number of pairs.", "\n \n")
cat("McNemar's Chi^2 Statistic (corrected for continuity) = ", round(object$XMc, digits=object$digits), " which has a p-value of: ", round(object$XMc.p.value, digits=object$digits), "\n",sep="")
cat("Note: The p.value for McNemar's Test corresponds to the hypothesis test: H0: OR = 1 vs. HA: OR != 1", "\n", sep="")

cat("McNemar's Odds Ratio (b/c): ", round(object$ORMc, digits=object$digits), "\n", sep="") 
cat(100*(1-object$alpha), "% Confidence Limits for the OR are: [", round(object$ORMc.CIL,digits=object$digits), ", ", round(object$ORMc.CIU,digits=object$digits),"]", "\n", sep="");

cat("The risk difference is: ", round(object$rd, digits=object$digits), "\n", sep="") 
cat(100*(1-object$alpha), "% Confidence Limits for the rd are: [", round(object$rd.CIL,digits=object$digits), ", ", round(object$rd.CIU,digits=object$digits),"]", "\n", sep="");
}