epiKappa <- function(C, alpha=0.05, k0=0.4, digits=3)
{
#Error Checking
if (nrow(C) != ncol(C)) 
        stop("Sorry, this function Requires a Square Matrix of Counts of Agreement or Proportions of Agreement")

if ((alpha >= 1) || (alpha <= 0))
        stop("Sorry, the alpha must lie within (0,1)")

#Append Row/Col classifications
if ( (is.null(rownames(C))) && (is.null(colnames(C))) )
{
cN <- NULL; 
rN <- NULL;
for (i in 1:nrow(C))
{
cN <- append(cN, paste("Rater I: Type ", i, sep=""));
rN <- append(rN, paste("Rater II: Type ", i, sep=""));
}
colnames(C) <- cN
rownames(C) <- rN
}

#Initialize Objects
r <- NULL;
r$k0 <- k0; r$digits <- digits;
r$alpha <- alpha;
r$Data <- C;

T1 <- sum(C);
X <- C/T1;

r$p0 <- sum(diag(X));
pe <- rep(NA, nrow(X));

for (i in 1:nrow(X))
{
pe[i] <- sum(X[,i])*sum(X[i,]);
}

r$pe <- sum(pe)

#Calculate Kappa using pe
r$kappa <- (r$p0 - r$pe)/(1 - r$pe);

#Strings for Fleiss Agreement Statement
if (r$kappa <= 0.4)
{
r$Fleiss <- "poor";
}

if ( (0.4 < r$kappa) && (r$kappa <= 0.75) )
{
r$Fleiss <- "fair to good";
}

if ( (0.75 < r$kappa) && (r$kappa <= 1.0) )
{
r$Fleiss <- "excellent";
}

#Placeholders to calculate variance using standard formulae.
AH <- rep(NA, nrow(X)); 
BH <- matrix(0, nrow(X), nrow(X)); 
AC <- rep(NA, nrow(X)); 
BC <- matrix(0, nrow(X), nrow(X)); 

for (i in 1:nrow(X))
{
AH[i] <- X[i,i]*(( 1 - (sum(X[i,]) + sum(X[,i])) *(1-r$k0)))^2;
AC[i] <- X[i,i]*(( 1 - (sum(X[i,]) + sum(X[,i])) *(1-r$kappa)))^2;
}

for (i in 1:nrow(X))
{
for (j in 1:nrow(X))
{
if (i != j)
{
BH[i,j] <- (1 - r$k0)^2 * X[i,j]*(sum(X[j,]) + sum(X[,i]));
BC[i,j] <- (1 - r$kappa)^2 * X[i,j]*(sum(X[j,]) + sum(X[,i]));
}
}
}

AH <- sum(AH); AC <- sum(AC);
BH <- sum(BH); BC <- sum(BC);
CH <- (r$k0 -r$pe*(1-r$k0))^2; 
CC <- (r$kappa -r$pe*(1-r$kappa))^2;

#Standard errors for hypothesis tests and CI's
r$seH <- sqrt( (AH + BH -CH)/T1 )/(1-r$pe);
r$seC <- sqrt( (AC + BC -CC)/T1 )/(1-r$pe);

#Z Test for kappa
r$Z <- (r$kappa - r$k0)/r$seH;
r$p.value <- 1 - pnorm(r$Z);

#CI's for kappa.
r$CIL <- r$kappa - qnorm(1 - alpha/2) * r$seC;
r$CIU <- r$kappa + qnorm(1 - alpha/2)* r$seC;

class(r) <- "epiKappa";
return(r);
}

#summary method
summary.epiKappa <- function(object, ...)
{
cat("Kappa Analysis of Agreement", "\n \n")

print(object$Data);
cat("\n");

cat("Cohen's Kappa is: ", round(object$kappa, digits = object$digits), "\n")
cat("According to Fleiss (1981), the point estimate of kappa suggests ", object$Fleiss, " agreement.", "\n \n",sep="");

cat(100*(1-object$alpha), "% Confidence Limits for the true Kappa Statistic are: [", round(object$CIL,digits=object$digits), ", ", round(object$CIU,digits=object$digits),"]", "\n \n", sep="");
cat("Z Test for H0: kappa = ", object$k0, " vs. HA: kappa >= ", object$k0, " is ", round(object$Z,digits=object$digits), " with a p.value of ", round(object$p.value, digits=object$digits), "\n \n", sep="");
cat("The associated standard error under H0 is: ", round(object$seH, digits=object$digits), "\n", sep="")
}

#Print method
print.epiKappa <- function(x, ...)
{
cat("Kappa Analysis of Agreement", "\n \n")

cat("Cohen's Kappa is: ", round(x$kappa, digits = x$digits), "\n")
cat("According to Fleiss (1981), the point estimate of kappa suggests ", x$Fleiss, " agreement. \n", sep="");
}