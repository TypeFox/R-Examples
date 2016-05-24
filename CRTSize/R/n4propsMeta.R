n4propsMeta <- function(data, k, ICC, ICCDistn="unif", lower=0, upper=0.25, varRed=FALSE, m, sdm, pC, sdpC, iter=1000, alpha=0.05, RR=TRUE)
{
if (!is.matrix(data))
	stop("Sorry data must be a matrix of RR/OR, 95 % Lower and Upper Limits from Previous Studies.")

if (! ( (ICCDistn == "fixed") | (ICCDistn == "unif") | (ICCDistn == "normal") | (ICCDistn == "smooth") ) )
	stop("Sorry, the ICC Distribution must be one of: fixed, unif, normal or smooth.")

if ( (ICCDistn == "fixed") && (length(ICC) > 1) )
	stop("Sorry you can only provide a single ICC value with the fixed distribution option.")

if ((alpha >= 1) || (alpha <= 0))
        stop("Sorry, the alpha must lie within (0,1).")

if ((sdm < 0) || (sdpC < 0))
        stop("Sorry, the standard deviations of m must be non-negative.")

if ( (pC <= 0) || (pC >= 1) )
        stop("Sorry, the Control Rate must lie within (0,1).")

for (i in 1:length(ICC))
{
if (ICC[i] <= 0)
        stop("Sorry, the ICC must lie within (0,1).")
}

if (m <=1) 
        stop("Sorry, the (average) cluster size, m, should be greater than one.")

for (i in 1:length(k))
{
if (k[i] <= 1)
        stop("Sorry, the values of k must be greater than 1")
}


X <- NULL;

X$data <- data; X$k <- k; X$ICC <- ICC;
X$m <- m; X$sdm <- sdm; X$pC <- pC;
X$sdpC <- sdpC; X$iter <- iter; X$ICCDistn <- ICCDistn;
X$lower <- lower; X$upper <- upper;
X$alpha <- alpha; X$varRed <- varRed;

original <- fixedMetaAnalRROR(data, alpha=X$alpha);

X$newMean <- original$thetaF;
X$newVar <- original$Var;

X$thetaNew <- rnorm(1, X$newMean, sqrt(X$newVar))
X$lF <- original$lF;
X$uF <- original$uF

#Obtain pT from the new simulated value and pC...

X$Power <- rep(0, length(k));

if (varRed)
{
X$varianceReduction <- rep(0, length(k));
}

for (a in 1:length(k))
{
kT0 <- k[a];
kC0 <- k[a];

if (ICCDistn == "unif")
{
ICCT0 <- runif(iter, lower, upper)
}

if (ICCDistn == "fixed")
{
ICCT0 <- rep(ICC, iter)
}

if (ICCDistn == "normal")
{
ICCT0 <- abs( rnorm(iter, 0, sd(ICC) ) );
}

if (ICCDistn == "smooth")
{
dens <- density(ICC, n=2^16, from=lower, to=upper)
ICCT0 <- sample(dens$x, size=iter, prob=dens$y)
}

Reject <- rep(NA, iter);

if (varRed)
{
varReductionIter <- rep(NA,iter);
}

for (i in 1:iter)
{
pC0 <- rnorm(1, pC, sdpC);

pT0 <- exp(X$thetaNew + log(pC0));

w <- .oneCRTBinary(pT=pT0, pC=pC0, kC=kC0, kT=kT0, mTmean=m, mTsd=sdm, mCmean=m, mCsd=sdm, ICCT=ICCT0[i], ICCC=ICCT0[i])
x <- .summarizeTrialRROR(ResultsTreat=w$ResultsTreat, ResultsControl=w$ResultsControl, RR=RR)
y <- .makeCIRROR(logRROR=x$logRROR, varlogRROR=x$VarLogRROR)
z <- fixedMetaAnalRROR(data=rbind(data, y), alpha=alpha);
Reject[i] <- z$Sig;

if (varRed)
{
varReductionIter[i] <- z$Var;
}

}

X$Power[a] <- sum(Reject, na.rm=TRUE)/iter;

if (varRed)
{
X$varianceReduction[a] <- mean(varReductionIter, na.rm=TRUE)/X$newVar;
}
}

names(X$Power) <- k

if (varRed)
{
names(X$varianceReduction) <- k
}

class(X) <- "n4propsMeta";
return(X);

}


#Print Method
print.n4propsMeta <- function(x, ...)
{
cat("The Approximate Power of the Updated Meta-Analysis is: (Clusters per Group) \n");
print(x$Power);

if (x$varRed)
{
cat("The Approximate Proportion of Variance Reduction is: (Clusters per Group) \n");
print(1 - x$varianceReduction);
}

}

#Summary method
summary.n4propsMeta <- function(object, ...)
{
cat("Sample Size Calculation for Binary Outcomes Based on Updated Meta-Analysis", "\n \n", sep="")
cat("The original fixed effects Relative Risk/Odds Ratio is ", exp(object$newMean), "\n", sep="");
cat("With ", (1 - object$alpha)*100,  "% Confidence Limits: (", exp(object$lF), ", ", exp(object$uF), ") \n \n",sep="");

cat("The Approximate Power of the Updated Meta-Analysis is: (Clusters per Group) \n", sep="");
print(object$Power);

if (object$varRed)
{
cat("The Approximate Proportion of Variance Reduction is: (Clusters per Group) \n");
print(1 - object$varianceReduction);
}
cat("\n", "Assuming:", "\n", sep="")
cat("Proportion with Outcome in Control Group: ", object$pC, " with standard deviation: ", object$sdpC,  "\n", sep="");
cat("Mean Cluster Size: ", object$m, " with standard deviation: ", object$sdm, "\n", sep="");
cat("ICC =", object$ICC, "\n");
cat("ICC Distribution", object$ICCDistn, "\n");
cat("Clusters =", object$k, "\n");
cat("Iterations =", object$iter, "\n");
}

#############################################
#A couple of basic helper functions;

#Takes a trial, from oneCRT function; #Generates RR or OR and variances;

.summarizeTrialRROR <- function(ResultsTreat, ResultsControl, RR=TRUE)
{
Summary <- NULL;

if (RR)
{
Summary$RROR <- ResultsTreat[1]/ResultsControl[1];
Summary$logRROR <- log(Summary$RROR);
Summary$VarLogRROR <- ( ((1 - ResultsTreat[1])*ResultsTreat[3])/(ResultsTreat[2]*ResultsTreat[1]) + ((1 - ResultsControl[1])*ResultsControl[3])/(ResultsControl[2]*ResultsControl[1]) )
}

else
{
Summary$RROR <- (ResultsTreat[1]/(1 - ResultsTreat[1]))/ (ResultsControl[1]/(1 - ResultsControl[1]));
Summary$logRROR <- log(Summary$RROR);
Summary$VarLogRROR <- ResultsTreat[3]/(ResultsTreat[2]*ResultsTreat[1]*(1 - ResultsTreat[1])) + ResultsControl[3]/(ResultsControl[2]*ResultsControl[1]*(1 - ResultsControl[1]))
}
return(Summary);
}

###############################

#Returns Confidence Interval for either OR or RR
.makeCIRROR <- function(logRROR, varlogRROR)
{
X <- c(exp(logRROR), exp(logRROR - 1.96*sqrt(varlogRROR)), exp(logRROR + 1.96*sqrt(varlogRROR)))
return(X);
}

###################

#Internal method for generating clustered binary data according to Lunn and Davies;


.oneCRTBinary <- function(pC, pT, kC, kT, mTmean, mTsd, mCmean, mCsd, ICCT, ICCC)
{
X <- NULL;

#Treatment Loop, generate

mT <- floor(rnorm(kT, mTmean, mTsd))

for (j in 1:kT)
{
if (mT[j] <= 10)
{
mT[j] <- 10;
}
}


dataT <- matrix(NA, nrow=max(mT), ncol=kT);

for (j in 1:kT)
{
Z <- rbinom(1, 1, pT);

for (i in 1:mT[j])
{

U <- rbinom(1, 1, sqrt(ICCT));
Y <- rbinom(1, 1, pT);
dataT[i,j] <- (1 - U)*Y + U*Z; 

}
}


#Control Loop...
 
mC <- floor(rnorm(kC, mCmean, mCsd))

for (j in 1:kC)
{
if (mC[j] <= 10)
{
mC[j] <- 10;
}
}

dataC <- matrix(NA, nrow=max(mC), ncol=kC);

for (j in 1:kC)
{
Z <- rbinom(1, 1, pC);

for (i in 1:mC[j])
{

U <- rbinom(1, 1, sqrt(ICCC));
Y <- rbinom(1, 1, pC);
dataC[i,j] <- (1 - U)*Y + U*Z; 

}
}

#Total number in treatment and control groups;
MTreat <- nrow(dataT)*ncol(dataT) - sum(is.na(dataT))
MControl <- nrow(dataC)*ncol(dataC) - sum(is.na(dataC))

#PhatT is the treatment rate;
PhatT <- sum(dataT, na.rm=TRUE)/MTreat

if ((PhatT == 0) || is.na(PhatT) || is.infinite(PhatT) )
{
PhatT <- 1/MTreat;
}

#PhatC is control rate;
PhatC <- sum(dataC, na.rm=TRUE)/MControl


if ((PhatC == 0) || is.na(PhatC) || is.infinite(PhatC) )
{
PhatC <- 1/MControl;
}

#Average cluster size;
mbarT <- sum(mT^2)/sum(mT);
mbarC <- sum(mC^2)/sum(mC);

#ICC Calculations

MSC <- 0; MSW <- 0;

for (j in 1:kT)
{
MSC <- MSC + mT[j]*( sum(dataT[,j], na.rm=TRUE)/mT[j] -  PhatT)^2 / (kT + kC - 2);
MSW <- MSW + mT[j]*(sum(dataT[,j], na.rm=TRUE)/mT[j])*(1 - (sum(dataT[,j], na.rm=TRUE)/mT[j]))/ (MTreat + MControl -(kT + kC));
}

for (j in 1:kC)
{
MSC <- MSC + mC[j]*(sum(dataC[,j], na.rm=TRUE)/mC[j] -  PhatC)^2 / (kT + kC - 2);
MSW <- MSW + mC[j]*(sum(dataC[,j], na.rm=TRUE)/mC[j])*(1 - (sum(dataC[,j], na.rm=TRUE)/mC[j]))/ (MTreat + MControl -(kT + kC));
}

#Assume a common ICC between treatment and control groups;
m0 <- ( (MTreat + MControl) - (mbarT + mbarC) ) / ( (kT + kC) - 2);
ICC <- max((MSC - MSW)/(MSC + (m0 - 1)*MSW), 0);

#Inflation factors;
CT <- (1 + (mbarT - 1)*ICC);
CC <- (1 + (mbarC - 1)*ICC);

#Summary of what is typically required, including rates, total number of subjects (group), clusters and ICC;
X$ResultsTreat <- c(PhatT, MTreat, CT, ICC, kT);
X$ResultsControl <- c(PhatC, MControl, CC, ICC, kC);

return(X);
}