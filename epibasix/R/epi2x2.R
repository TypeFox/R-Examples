epi2x2 <- function(X, alpha=0.05, digits=3)
{
#Error Checking
if (nrow(X) != ncol(X)) 
        stop("Sorry, 'X' must be a 2x2 Matrix...")

if ((alpha >= 1) || (alpha <= 0))
        stop("Sorry, the alpha must lie within (0,1)")

#Assign Row and Column Labelss
if ( (is.null(rownames(X))) && (is.null(colnames(X))) )
{
rownames(X) <- c("Risk Present", "Risk Absent")
colnames(X) <- c("Disease Present (Cases)", "Disease Absent (Controls)");
}

#Initialize Parameters
r <- NULL;
r$digits <- digits;
r$X <- X;
r$alpha <- alpha;

#For notation, r1 = sum of first row, etc.
a <- X[1,1]; b <- X[1,2]; c <- X[2,1]; d <- X[2,2];
r1 <- a + b; r2 <- c + d; c1 <- a + c; c2 <- b+d;
T1 <- r1 + r2;

#Chi-Squared Test (With Continuity Correction)
r$Sy <- T1*((abs(a*d - b*c) - T1/2)^2)/(c1*c2*r1*r2);
r$Sy.p.value <- 1- pchisq(r$Sy,1);  

#Fisher's Exact
r$Fisher.p.value <- fisher.test(X)$p.value;

#Odds Ratio and Confidence Limits
r$OR <- (a*d)/(b*c)

r$OR.CIL <- exp(log(r$OR) - qnorm(1 - alpha/2) * sqrt(1/a + 1/b + 1/c + 1/d));
r$OR.CIU <- exp(log(r$OR) + qnorm(1 - alpha/2) * sqrt(1/a + 1/b + 1/c + 1/d));

#Cohort Study Info
r$p1Co <- a/r1;  r$p2Co <- c/r2;

r$rdCo <- r$p1Co - r$p2Co;
r$rdCo.CIL <- r$rdCo - (1/r1 + 1/r2)/2 - qnorm(1 - alpha/2) * sqrt(r$p1Co*(1-r$p1Co)/r1 + r$p2Co*(1-r$p2Co)/r2);
r$rdCo.CIU <- r$rdCo + (1/r1 + 1/r2)/2 + qnorm(1 - alpha/2) * sqrt(r$p1Co*(1-r$p1Co)/r1 + r$p2Co*(1-r$p2Co)/r2);

#Relative Risk and CI
r$RR <- r$p1Co/r$p2Co;
r$RR.CIL <- exp( log(r$RR) - qnorm(1 - alpha/2) * sqrt( b/(a*(a+b)) + d/(c*(c+d)) ));
r$RR.CIU <- exp( log(r$RR) + qnorm(1 - alpha/2) * sqrt( b/(a*(a+b)) + d/(c*(c+d)) ));

#Case Control Study risk difference
r$p1CC <- a/c1; r$p2CC <- b/c2;
r$rdCC <- r$p1CC - r$p2CC;
r$rdCC.CIL <- r$rdCC - (1/c1 + 1/c2)/2 - qnorm(1 - alpha/2) * sqrt(r$p1CC*(1-r$p1CC)/r1 + r$p2CC*(1-r$p2CC)/r2);
r$rdCC.CIU <- r$rdCC + (1/c1 + 1/c2)/2 + qnorm(1 - alpha/2) * sqrt(r$p1CC*(1-r$p1CC)/r1 + r$p2CC*(1-r$p2CC)/r2);

class(r) <- "epi2x2";
return(r);
}

#Print Method
print.epi2x2 <- function(x, ...)
{
cat("Epidemiological 2x2 Table Analysis", "\n \n");
cat("Input Matrix: \n")
print(x$X)
cat("\n", "Pearson Chi-Squared Statistic (Includes Yates' Continuity Correction): ", round(x$Sy, digits=x$digits), "\n", sep="")
cat("Associated p.value for H0: There is an association between exposure and outcome vs. HA: No association : ",  round(x$Sy.p.value, digits=x$digits), "\n", sep="")
cat("p.value using Fisher's Exact Test: ",  round(x$Fisher.p.value, digits=x$digits), "\n \n", sep="");

cat("Estimate of Odds Ratio: ",  round(x$OR, digits=x$digits), "\n", sep="")
cat(100*(1-x$alpha), "% Confidence Limits for true Odds Ratio are: [", round(x$OR.CIL,digits=x$digits), ", ", round(x$OR.CIU,digits=x$digits),"]", "\n", sep="");
}

#Summary Method
summary.epi2x2 <- function(object, ...)
{
cat("Epidemiological 2x2 Table Analysis", "\n \n");
cat("Input Matrix: \n")
print(object$X)
cat("\n", "Pearson Chi-Squared Statistic (Includes Yates' Continuity Correction): ", round(object$Sy, digits=object$digits), "\n", sep="")
cat("Associated p.value for H0: There is an association between exposure and outcome vs. HA: No association : ",  round(object$Sy.p.value, digits=object$digits), "\n", sep="")
cat("p.value using Fisher's Exact Test (1 DF) : ",  round(object$Fisher.p.value, digits=object$digits), "\n \n", sep="");

cat("Estimate of Odds Ratio: ",  round(object$OR, digits=object$digits), "\n", sep="")
cat(100*(1-object$alpha), "% Confidence Limits for true Odds Ratio are: [", round(object$OR.CIL,digits=object$digits), ", ", round(object$OR.CIU,digits=object$digits),"]", "\n \n", sep="");

cat("Estimate of Relative Risk (Cohort, Col1): ",  round(object$RR, digits=object$digits), "\n", sep="")
cat(100*(1-object$alpha), "% Confidence Limits for true Relative Risk are: [", round(object$RR.CIL,digits=object$digits), ", ", round(object$RR.CIU,digits=object$digits),"]", "\n \n", sep="");

cat("Estimate of Risk Difference (p1 - p2) in Cohort Studies: ",  round(object$rdCo, digits=object$digits), "\n", sep="")
cat(100*(1-object$alpha), "% Confidence Limits for Risk Difference: [", round(object$rdCo.CIL,digits=object$digits), ", ", round(object$rdCo.CIU,digits=object$digits),"]", "\n \n", sep="");

cat("Estimate of Risk Difference (p1 - p2) in Case Control Studies: ",  round(object$rdCC, digits=object$digits), "\n", sep="")
cat(100*(1-object$alpha), "% Confidence Limits for Risk Difference: [", round(object$rdCC.CIL,digits=object$digits), ", ", round(object$rdCC.CIU,digits=object$digits),"]", "\n \n", sep="");

cat("Note: Above Confidence Intervals employ a continuity correction.", "\n")
}