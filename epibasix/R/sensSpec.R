sensSpec <- function(X, alpha=0.05, CL=TRUE, digits=3)
{
#Error Checking
if ((nrow(X) != ncol(X)) || (nrow(X) != 2))
        stop("Sorry, 'X' must be a 2x2 Matrix...")

if ((alpha >= 1) || (alpha <= 0))
        stop("Sorry, the alpha must lie within (0,1)")

#Assign row and column names
if ( (is.null(rownames(X))) && (is.null(colnames(X))) )
{
rownames(X) <- c("Reported A", "Reported B")
colnames(X) <- c("Gold Standard A", "Gold Standard B");
}

#Object Initialization
results <- NULL;
results$X <- X; results$digits <- digits;
results$alpha <- alpha;
results$CL <- CL;

#For notation, r1 = sum of first row, c1 = sum of first col. etc.
a <- X[1,1]; b <- X[1,2]; c <- X[2,1]; d <- X[2,2];
r1 <- a + b; r2 <- c + d; c1 <- a + c; c2 <- b+d;
T1 <- r1 + r2;

#Calculate point estimates...
results$sens <- a/c1; results$spec <- d/c2;

results$YoudenJ <- results$sens + results$spec - 1;

results$PA <- (a +d)/(T1);

#Compute Confidence Limits if required
if (CL)
{
results$sens.s <- sqrt(((results$sens*(1 - results$sens))/c1))
results$sens.CIL <- results$sens - qnorm(1 - alpha/2)*results$sens.s
results$sens.CIU <- results$sens + qnorm(1 - alpha/2)*results$sens.s

results$spec.s <- sqrt(((results$spec*(1 - results$spec))/c2))
results$spec.CIL <- results$spec - qnorm(1 - alpha/2)*results$spec.s
results$spec.CIU <- results$spec + qnorm(1 - alpha/2)*results$spec.s

results$PA.s <- sqrt(((results$PA*(1 - results$PA))/T1))
results$PA.CIL <- results$PA - qnorm(1 - alpha/2)*results$PA.s
results$PA.CIU <- results$PA + qnorm(1 - alpha/2)*results$PA.s

results$YoudenJ.s <- sqrt( ((results$sens*(1 - results$sens))/c1) + ((results$spec*(1 - results$spec))/c2) );
results$YoudenJ.CIL <- results$YoudenJ - qnorm(1 - alpha/2)*results$YoudenJ.s
results$YoudenJ.CIU <- results$YoudenJ + qnorm(1 - alpha/2)*results$YoudenJ.s
}

class(results) <- "sensSpec";

return(results);

}

#Print Method
print.sensSpec <- function(x, ...)
{
cat("\n", "Simple Sensitivity and Specitivity Output", "\n \n");
cat("Input Matrix: \n")
print(x$X)
cat("\n", "The sample of sensitivity is: ", 100* round(x$sens,digits=x$digits), "% \n", sep="") 
cat("\n", "The sample of specificity is: ", 100* round(x$spec,digits=x$digits), "% \n", sep="") 
}

#Summary Method
summary.sensSpec <- function(object, ...)
{
cat("\n", "Detailed Sensitivity and Specitivity Output", "\n \n");

cat("Input Matrix: \n")
print(object$X)
cat("\n", "The sample sensitivity is: ", 100* round(object$sens,digits=object$digits), "% \n", sep="") 
if (object$CL){
cat(100*(1-object$alpha), "% Confidence Limits for true sensitivity are: [", 100* round(object$sens.CIL,digits=object$digits), ", ", 100* round(object$sens.CIU,digits=object$digits),"]", "\n", sep="");
}
cat("\n", "The sample of specificity is: ", 100* round(object$spec,digits=object$digits), "% \n", sep="") 
if (object$CL){
cat(100*(1-object$alpha), "% Confidence Limits for true specificity are: [",100*  round(object$spec.CIL,digits=object$digits), ", ", 100* round(object$spec.CIU,digits=object$digits),"]", "\n", sep="");
}


cat("\n", "The sample value of Youden's J is: ", 100*round(object$YoudenJ,digits=object$digits), "\n", sep="") 
if (object$CL){
cat(100*(1-object$alpha), "% Confidence Limits for Youden's J are: [", 100*round(object$YoudenJ.CIL,digits=object$digits), ", ", 100*round(object$YoudenJ.CIU,digits=object$digits),"]", "\n", sep="");
}
cat("\n","Sample value for Percent Agreement (PA) is: ",  100*round(object$PA,digits=object$digits), "% \n", sep="")
if (object$CL){
cat(100*(1-object$alpha), "% Confidence Limits for PA are: [", 100*round(object$PA.CIL,digits=object$digits), ", ", 100*round(object$PA.CIU,digits=object$digits),"]", "\n \n", sep="");
}
}