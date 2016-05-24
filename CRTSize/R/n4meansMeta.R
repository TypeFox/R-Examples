n4meansMeta <- function(data, k, ICC, ICCDistn="unif", lower=0, upper=0.25, varRed=FALSE, m, sdm, meanC, sdC, sdT=sdC, iter=1000, alpha=0.05)
{
if (!is.matrix(data))
	stop("Sorry data must be a matrix of Mean Difference, 95 % Lower and Upper Limits from Previous Studies")

if (! ( (ICCDistn == "fixed") | (ICCDistn == "unif") | (ICCDistn == "normal") | (ICCDistn == "smooth") ) )
	stop("Sorry, the ICC Distribution must be one of: fixed, unif, normal or smooth.")

if ( (ICCDistn == "fixed") && (length(ICC) > 1) )
	stop("Sorry you can only provide a single ICC value with the fixed distribution option.")

if ((alpha >= 1) || (alpha <= 0))
        stop("Sorry, the alpha must lie within (0,1)")

if ((sdm < 0) || (sdC < 0))
        stop("Sorry, variances must be non-negative.")


for (i in 1:length(ICC))
{
if (ICC[i] <= 0)
        stop("Sorry, the ICC must lie within (0,1)")
}

if (m <=1) 
        stop("Sorry, the (average) cluster size, m, should be greater than one...")

for (i in 1:length(k))
{
if (k[i] <= 1)
        stop("Sorry, the values of k must be greater than 1")
}


X <- NULL;

X$data <- data; X$k <- k; X$ICC <- ICC;
X$m <- m; X$sdm <- sdm; X$meanC <- meanC;
X$sdC <- sdC; X$iter <- iter; X$varRed <- varRed;
X$alpha <- alpha;

original <- fixedMetaAnalMD(data, alpha=X$alpha);

X$newMean <- original$thetaF;
X$newVar <- original$Var;

X$thetaNew <- rnorm(1, X$newMean, sqrt(X$newVar))
X$lF <- original$lF;
X$uF <- original$uF

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
meanC0 <- rnorm(1, meanC, sdC);

meanT0 <- X$thetaNew + meanC0;

w <- .oneCRTCTS(meanT=meanT0, meanC=meanC0, sdT=sdT, sdC=sdC, kC=kC0, kT=kT0, mTmean=m, mTsd=sdm, mCmean=m, mCsd=sdm, ICCT=ICCT0[i], ICCC=ICCT0[i])
x <- .summarizeTrialCTS(ResultsTreat=w$ResultsTreat, ResultsControl=w$ResultsControl)
y <- .makeCI(Delta=x$meanDiff, varDelta=x$varMeanDiff)
z <- fixedMetaAnalMD(data=rbind(data, y), alpha=alpha);
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

class(X) <- "n4meansMeta";
return(X);

}


#Print Method
print.n4meansMeta <- function(x, ...)
{
cat("The Approximate Power of the Updated Meta-Analysis is: (Clusters per Group) \n");
print(x$Power);
if (x$varRed)
{
cat("The Approximate Proportion of Variance Reduction is: (Clusters per Group) \n");
print(1 - x$varianceReduction);
}
}


#Print method
summary.n4meansMeta <- function(object, ...)
{
cat("Sample Size Calculation for Binary Outcomes Based on Updated Meta-Analysis", "\n \n", sep="")
cat("The original Fixed Effects Risk Difference is ", object$newMean, "\n", sep="");
cat("With ", (1 - object$alpha)*100,  "% Confidence Limits: (", object$lF, ", ", object$uF, ") \n \n",sep="");

cat("The Approximate Power of the Updated Meta-Analysis is: (Clusters per Group) \n", sep="");
print(object$Power);

if (object$varRed)
{
cat("The Approximate Proportion of Variance Reduction is: (Clusters per Group) \n");
print(1 - object$varianceReduction);
}

cat("\n", "Assuming:", "\n", sep="")
cat("Control Group Mean: ", object$meanC, " with standard deviation: ", object$sdC,  "\n", sep="");
cat("Mean Cluster Size: ", object$m, " with standard deviation: ", object$sdm, "\n", sep="");
cat("ICC =", object$ICC, "\n");
cat("ICC Distribution", object$ICCDistn, "\n");
cat("Clusters =", object$k, "\n");
cat("Iterations =", object$iter, "\n");
}

#################################################
#A couple of basic helper functions;

#Summarize: Requires ResultsTreat and ResultsControl to give information about the mean diff and its variance
#X$ResultsTreat <- c(meanTreat, MTreat, CT, ICC, MSC, MSW, VarTreat,kT);
#X$ResultsControl <- c(meanControl, MControl, CC, ICC, MSC, MSW, VarControl, kC);

.summarizeTrialCTS <- function(ResultsTreat, ResultsControl, ...)
{
Summary <- NULL;

Summary$meanDiff <- ResultsTreat[1] - ResultsControl[1];
Summary$varMeanDiff <- (ResultsTreat[5] + ResultsTreat[6]);

return(Summary);
}


###############################

#Returns Confidence Interval
.makeCI <- function(Delta, varDelta)
{
X <- c(Delta, (Delta - 1.96*sqrt(varDelta)), (Delta + 1.96*sqrt(varDelta)))
return(X);
}

###################

#Simple function to count number of 2500's, as very large clusters
#cause the system to run out of memory.

.number2500s <- function(n)
{
x <- floor(n/2500);
return(x);
}


#####################
#Internal Method for generation of Clustered CTS Data (Normally Distributed)

.oneCRTCTS <- function(meanT, sdT, meanC, sdC, kT, kC, mTmean, mTsd, mCmean, mCsd, ICCT, ICCC)
{
X <- NULL;

####Treatment Loop, generate CTS data in blocks..

mT <- floor(rnorm(kT, mTmean, mTsd))

dataT <- matrix(NA, nrow=max(mT), ncol=kT);

for (j in 1:kT)
{

a <- .number2500s(mT[j]);

if (a >= 1)
{

for (k in 1:a)
{
Sigma <- matrix(ICCT*sdT^2, nrow=2500, ncol=2500);
diag(Sigma) <- c(rep(sdT^2, 2500));
dataT[((k-1)*2500 +1):(2500*k),j] <- round(.mvrnorm(n=1, mu=rep(meanT, 2500), Sigma=Sigma), digits=5)
}

if ( (mT[j] - 2500*a) != 0)
{
lastones <- mT[j] - 2500*a;
Sigma <- matrix(ICCT*sdT^2, nrow=lastones, ncol=lastones);
diag(Sigma) <- c(rep(sdT^2, lastones));
dataT[(2500*k+1):mT[j],j] <- round(.mvrnorm(n=1, mu=rep(meanT, lastones), Sigma=Sigma), digits=5)
}

}

if (a < 1)
{
lastones <- mT[j] - 2500*a;
Sigma <- matrix(ICCT*sdT^2, nrow=lastones, ncol=lastones);
diag(Sigma) <- c(rep(sdT^2, lastones));
dataT[1:mT[j],j] <- round(.mvrnorm(n=1, mu=rep(meanT, lastones), Sigma=Sigma), digits=5)
}

}

####Control Loop

mC <- floor(rnorm(kC, mCmean, mCsd))

dataC <- matrix(NA, nrow=max(mC), ncol=kC);

j <- 1;
k <- 1;

for (j in 1:kC)
{

a <- .number2500s(mC[j]);

if (a >= 1)
{
for (k in 1:a)
{
Sigma <- matrix(ICCC*sdC^2, nrow=2500, ncol=2500);
diag(Sigma) <- c(rep(sdC^2, 2500));
dataC[((k-1)*2500 +1):(2500*k),j] <- round(.mvrnorm(n=1, mu=rep(meanC, 2500), Sigma=Sigma), digits=5)
}

if ( (mC[j] - 2500*a) != 0)
{
lastones <- mC[j] - 2500*a;
Sigma <- matrix(ICCC*sdC^2, nrow=lastones, ncol=lastones);
diag(Sigma) <- c(rep(sdC^2, lastones));
dataC[(2500*k+1):mC[j],j] <- round(.mvrnorm(n=1, mu=rep(meanC, lastones), Sigma=Sigma), digits=5)
}

}
if (a < 1)
{
lastones <- mC[j] - 2500*a;
Sigma <- matrix(ICCC*sdC^2, nrow=lastones, ncol=lastones);
diag(Sigma) <- c(rep(sdC^2, lastones));
dataC[1:mC[j],j] <- round(.mvrnorm(n=1, mu=rep(meanC, lastones), Sigma=Sigma), digits=5)
}

}

meanTreat <- mean(dataT, na.rm=TRUE);
meanControl <- mean(dataC, na.rm=TRUE);


MTreat <- nrow(dataT)*ncol(dataT) - sum(is.na(dataT))
MControl <- nrow(dataC)*ncol(dataC) - sum(is.na(dataC))


mbarT <- sum(mT^2)/sum(mT);
mbarC <- sum(mC^2)/sum(mC);

m0 <- ( (MTreat + MControl) - (mbarT + mbarC) ) / ( (kT + kC) - 2);

#ICC Calculations

MSC <- 0; MSW <- 0;
VarTreat <- 0; VarControl <- 0;

for (j in 1:kT)
{
MSC <- MSC + mT[j]*( sum(dataT[,j], na.rm=TRUE)/mT[j] -  meanTreat)^2 / (kT + kC - 2);

for (i in 1:mT[j])
{
MSW <- MSW + ((dataT[i,j] - mean(dataT[,j], na.rm=TRUE))^2/ (MTreat + MControl -(kT + kC)));
VarTreat <- VarTreat + ((dataT[i,j] - meanTreat)^2/(MTreat - 1))
}

}

for (j in 1:kC)
{
MSC <- MSC + mC[j]*( sum(dataC[,j], na.rm=TRUE)/mC[j] -  meanControl)^2 / (kT + kC - 2);

for (i in 1:mC[j])
{
MSW <- MSW + ((dataC[i,j] - mean(dataC[,j], na.rm=TRUE))^2/ (MTreat + MControl - (kT + kC)));
VarControl <- VarControl + ((dataC[i,j] - meanControl)^2/(MControl - 1))
}

}


ICC <- max((MSC - MSW)/(MSC + (m0 - 1)*MSW), 0);

CT <- (1 + (mbarT - 1)*ICC);
CC <- (1 + (mbarC - 1)*ICC);


X$ResultsTreat <- c(meanTreat, MTreat, CT, ICC, MSC, MSW, VarTreat, kT);
X$ResultsControl <- c(meanControl, MControl, CC, ICC, MSC, MSW, VarControl, kC);

return(X);
}

.mvrnorm <- function (n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE) 
{
    p <- length(mu)
    if (!all(dim(Sigma) == c(p, p))) 
        stop("incompatible arguments")
    eS <- eigen(Sigma, symmetric = TRUE)
    ev <- eS$values
    if (!all(ev >= -tol * abs(ev[1]))) 
        stop("'Sigma' is not positive definite")
    X <- matrix(rnorm(p * n), n)
    if (empirical) {
        X <- scale(X, TRUE, FALSE)
        X <- X %*% svd(X, nu = 0)$v
        X <- scale(X, FALSE, TRUE)
    }
    X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% 
        t(X)
    nm <- names(mu)
    if (is.null(nm) && !is.null(dn <- dimnames(Sigma))) 
        nm <- dn[[1]]
    dimnames(X) <- list(nm, NULL)
    if (n == 1) 
        drop(X)
    else t(X)
}
