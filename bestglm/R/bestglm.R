bestglm <-
function (Xy, family=gaussian, IC = "BIC", t="default", CVArgs="default",  
    qLevel=0.99, TopModels=5, method="exhaustive", intercept=TRUE, weights=NULL, 
    nvmax="default", RequireFullEnumerationQ=FALSE, ...)
{
stopifnot(is.data.frame(Xy), ncol(Xy) > 1, 
    length(names(Xy)) == ncol(Xy), length(IC)>0)
if (any(is.na(match(IC, c("AIC","BIC","BICg","BICq", "CV", "LOOCV")))))
    stop(paste("IC =", IC, "invalid."))
if (any(is.na(match(method, c("exhaustive", "backward", "forward", "seqrep")))))
    stop(paste("method = ", method, "invalid."))
#
IncludeInterceptQ<-intercept
CVQ <- IC=="CV" || IC=="LOOCV"
if (!IncludeInterceptQ && CVQ)
    stop("Cross-validation not available with models without intercept term")
if (any(nvmax!="default") && CVQ)
    stop("Cross-validation not available with nvmax!=\"default\"")
if (is.character(family))
    stop("character string not allowed for argument 'family'")
binomialQ <- "binomial"==deparse(substitute(family))
gaussianQ <- "gaussian"==deparse(substitute(family))
#
y <- Xy[, ncol(Xy)]
n <- length(y)
#now check if binomial
LastColumnBinaryQ <- length(unique(y))==2
#
#in case of general binomial (non-logistic), adjust y to be the last 2 columns of Xy
SFQ <- FALSE
if (binomialQ && !LastColumnBinaryQ){
    SFQ <- TRUE
    p <- ncol(Xy)-2
    y<-SF<-as.matrix.data.frame(Xy[,c(p+1, p+2)])
    X<-Xy[,1:p]
    if (any(y<0))
        stop("Binomial nonlogistic-regression: S and F counts can not be <0")
    if (length(unique(apply(y, MARGIN=1, sum)))>1)
        stop("Binomial nonlogistic-regression: S + F must be constant")
}    
else {
    X <- Xy[, -ncol(Xy),drop=FALSE]
    p <- ncol(X)
    }
if (is.character(nvmax)){
    if (nvmax=="default")
        NVMAX <- p
    else
        stop(paste("error: either nvax = \"default\" or set to a number less than or equal to:",ncol(Xy[,-1])))
    }
else
    NVMAX <- min(nvmax, p)
#M - total number of subset models
    M <- 2^p-1
    if (TopModels > M)
        TopModels2 <- M
    else
        TopModels2 <- TopModels
    if (!is.na(match("CV",IC))||!is.na(match("LOOCV",IC)))
        TopModels2 <- 1
    if (TopModels2==1)
        BestModels <- NA
#At this point X and y are defined
#
#check that there is at least one input!
if (p==0)
    stop("dim(Xy) must be >=  2. Needs at least one covariate!")
#P is the number of terms in regression including the constant term
P <- p+1
#Xy0 needed for lm/glm
#Xy0 <- Xy
#names(Xy0)[P]<-"y"
#
#Determine settings for 'leaps' or Morgan-Tatar
#X dataframe with numeric or factor variables
if (any(unlist(lapply(X, is.character))))
    stop("input matrix must have numeric or factor variables only")
FactorsQ <- unlist(lapply(X, is.factor))
OrderedQ <- unlist(lapply(X, is.ordered))
NumDF <- rep(1, p) #initialize
#adjust NumDF for factor variables,
for (j in 1:ncol(X))
    if (FactorsQ[j])
        NumDF[j] <- length(levels(X[,j]))-1
TotalSize <- sum(NumDF) #total number of parameters (intercept not included)
CategoricalQ <- any(FactorsQ)
BinaryCategoricalQ <- CategoricalQ & all(NumDF==1)
if (any(NumDF>1) && CVQ)
    stop("Cross-validation not available when there are categorical variables with more than 2 levels!")
#if gaussian family, no factor with more than 2 levels, use 'leaps'
LEAPSQ <- gaussianQ && !any(NumDF>1) &&!RequireFullEnumerationQ
if (!LEAPSQ && !(IC=="BICq" && (t==0 || t==1))) {
    if (!gaussianQ)
        cat("Morgan-Tatar search since family is non-gaussian.",fill=TRUE)
    if (any(NumDF>1)&&gaussianQ)
        cat("Morgan-Tatar search since factors present with more than 2 levels.",fill=TRUE)
    if (any(NumDF>1)&&!gaussianQ)
        cat("Note: factors present with more than 2 levels.",fill=TRUE)
    if (RequireFullEnumerationQ&&gaussianQ&&all(NumDF==1))
        cat("Morgan-Tatar search RequireFullEnumerationQ=TRUE",fill=TRUE)
    }
#!LEAPSQ means Morgan-Tatar search. 
#In this case either 'lm' or #'glm' is used. 
#When FAMILY not "gaussian', 'glm' used. 
glmQ <- !(LEAPSQ || gaussianQ)
#Completed settings.
#
#Test if there is a column of 1's in X.
#Only do this test if:
#   1) X contains at least one column corresponding to a quantitative variable
#   2) IncludeInterceptQ (intercept is included already)
if (!all(FactorsQ)&& IncludeInterceptQ){
    X2<-X[,!FactorsQ,drop=FALSE] #removing factor variables
    p2<-ncol(X2)
    Column1Test<-apply(X2-matrix(rep(1,n*p2),ncol=p2), MARGIN=2, function(x) 0==sum(abs(x)))
    if (any(Column1Test))
      stop("Column of 1's can't be used! Intercept always included.")
}
#null deviance and df
if (glmQ)
    if (IncludeInterceptQ)
        ans <- glm(y ~ 1, family=family, weights=weights, ...)
    else 
        ans <- glm(y ~ -1, family=family, weights=weights, ...)
else
    if (IncludeInterceptQ)
        ans <- lm(y ~ 1, weights=weights,  ...)
    else
        ans <- lm(y ~ -1, weights=weights,  ...)
NullDeviance <- deviance(ans)
NullDF <- ans$df.residual
NullModel <- c(NullDeviance, NullDF)
names(NullModel) <- c("Deviance", "DF")
ModelReport <- list(NullModel=NullModel, LEAPSQ=LEAPSQ, glmQ=glmQ, gaussianQ=gaussianQ, NumDF=NumDF, 
    CategoricalQ=CategoricalQ, IncludeInterceptQ=IncludeInterceptQ)
#
#
#form ICLabel
if (!is.na(match("AIC",IC))) 
    ICLabel <- "AIC"
if (!is.na(match("BIC",IC))) 
    ICLabel <- "BIC"
if (!is.na(match("BICg",IC))) {
    if (is.character(t))
        g <- 1
    else
        g <- t
    ICLabel <- paste("BICg(g = ",g,")",sep="")
    }
if (!is.na(match("BICq",IC))) { 
    if (is.character(t))
        q <- 0.25
    else
        q <- t
    ICLabel <- paste("BICq(q = ",q,")",sep="")
    }
if (!is.na(match("CV",IC))) {
    if (is.list(CVArgs)){
        CVMethod<-CVArgs$Method
        K<-CVArgs$K
        REP<-CVArgs$REP
        stopifnot(!(is.null(CVMethod)|is.null(K)|is.null(REP)))
        }
    else if (CVArgs=="default"){
        CVMethod <- "d"
        K <- ceiling(n*(1-1/(log(n) - 1)))
        if (is.character(t))
            REP <- 1000
        else
            REP <- t
        }
    else
        stop("invalid CVArgs")
    if (any(is.na(match(CVMethod, c("HTF","DH","d")))))
        stop("invalid CVMethod")
    ICLabel <- switch(CVMethod,
        HTF= paste("CV(K = ", K, ", REP = ",REP, ")", sep=""),
        DH = paste("CVAdj(K = ", K,", REP = ",REP,  ")", sep=""),
        d = paste("CVd(d = ", K, ", REP = ",REP,  ")", sep="") )
}
if (IC=="CV" && CVMethod=="DH" && !gaussianQ)
    stop("DH cross-validation not available for non-Gaussian models.")
if (IC=="CV" && K < 2)
   stop("CV methods require K > 1!")
if (IC=="LOOCV" && !gaussianQ)
    stop("LOOCV cross-validation not available for non-Gaussian models.")
if (IC=="LOOCV") 
    ICLabel <- CVMethod <- "LOOCV"
#
#Ready for best regression subset search
#
#Now do we use exhaustive search
if (!LEAPSQ){
logLAllModels <- rep(-Inf,2^p)
AllModels <- rep(Inf,2^p) #initial IC value
#Special case -- we won't bother with full exhaustive enumeration
#  with BICq when t=0 (null model) and t=1 (full model)
if (glmQ)
    if (IncludeInterceptQ)
        ans <- glm(y ~ 1, family=family, weights=weights, ...)
    else 
        ans <- glm(y ~ -1, family=family, weights=weights, ...)
else
    if (IncludeInterceptQ)
        ans <- lm(y ~ 1, ...)
    else
        ans <- lm(y ~ -1, ...)
    if (IC=="BICq" && (t==0 || t==1)){
        if (t == 0) {
            if(glmQ) {
                    if (IncludeInterceptQ)
                        ans <- glm(y ~ 1, family=family, weights=weights, ...)
                else 
                    ans <- glm(y ~ -1, family=family, weights=weights, ...)
            }
            else {
                    if (IncludeInterceptQ)
                        ans <- lm(y ~ 1, weights=weights, ...)
                    else
                        ans <- lm(y ~ -1, weights=weights, ...)
            }
        }
        if (t == 1) {
            if(glmQ) {
                    if (IncludeInterceptQ)
                        ans <- glm(y ~ ., data=X, family=family, weights=weights, ...)
                else 
                    ans <- glm(y ~ -1+., data=X, family=family, weights=weights, ...)
            }
            else {
                    if (IncludeInterceptQ)
                        ans <- lm(y ~ ., data=X, weights=weights, ...)
                    else
                        ans <- lm(y ~ -1+., data=X, weights=weights, ...)
                }
            }
        out <- list(BestModel=ans, BestModels=NA, Bestq=NA, qTable=NA, Subsets=NA,    
            Title=ICLabel, ModelReport=ModelReport)
        class(out) <- "bestglm"
        cat("Note: in this special case with BICq with t =", t, "only fitted model is returned.",fill=TRUE)
        if (t ==0)
            cat("With t=0, null model is fitted.",fill=TRUE)
        else
            cat("With t=1, full model is fitted.",fill=TRUE)
        return(out)
    }
#Use glm and Morgan & Tatar method for exhaustive enumeration
# Best subsets by exhaustive search
# p>1
#Subsets has p rows indicating the presence/absence of each variable for best subset of size k=1,...,p. 
#There are p+1 columns including intercept.
#The null model, no variables apart from intercept, is included.
    logL<- rep(-Inf,p+1)
#Null model, no variable included
    kSize<- 0
    if (glmQ){
        if (IncludeInterceptQ)
            ans <- glm(y ~ 1, family=family, weights=weights, ...)
        else 
            ans <- glm(y ~ -1, family=family, weights=weights, ...)
            }
    else {
        if (IncludeInterceptQ)
            ans <- lm(y ~ 1, weights=weights, ...)
        else
            ans <- lm(y ~ -1, weights=weights, ...)
        }
    if (gaussianQ) #for consistency with 'leaps'
            L0 <- -(n/2)*log(sum(resid(ans)^2)/n)
    else
            L0 <- logLik(ans)
    Neg2LL <- -2*L0 
    penalty <- switch(IC,
                    AIC=2*kSize,
                    BIC=log(n)*kSize,
                    BICg=log(n)*kSize + 2*g*lchoose(TotalSize,kSize),
                    BICq=log(n)*kSize - 2*kSize*log(q/(1-q)),
                    CV = 0)
    BestModelIC <- Neg2LL + penalty
    AllModels[1] <- BestModelIC
    logLAllModels[1] <- L0
    BestModel<-ans
    if (IC=="BICq") {
            penalty <- Inf
            BestModelIC <- Inf
        }
    ICLabel3 <- paste("Sample Size: ",n,", Inputs: ",kSize, ", Penalty:",penalty, sep="")
#Search the best subsets over all non-null models
    for (i in 1:M){
        vars<-as.logical(rev(to.binary(i, p)))
        k<- sum(vars) #different from kSize!
        if (k > nvmax) 
            next
        Xi<-X[,vars,drop=FALSE]
        if (glmQ){
            if (IncludeInterceptQ)
                ans <- glm(y ~ ., data=Xi, family=family, weights=weights, ...)
            else 
                ans <- glm(y ~ -1+., family=family, weights=weights, ...)
            }
        else {
            if (IncludeInterceptQ)
                ans <- lm(y ~ ., data=Xi, weights=weights, ...)
            else
                ans <- lm(y ~ -1+., data=Xi, weights=weights, ...)
            }
        if (gaussianQ) #for consistency with 'leaps'
            L0 <- -(n/2)*log(sum(resid(ans)^2)/n)
        else
            L0 <- logLik(ans)
        kSize <- sum(NumDF[vars])
        Neg2LL <- -2*L0 
        penalty <- switch(IC,
                    AIC=2*kSize,
                    BIC=log(n)*kSize,
                    BICg=log(n)*kSize + 2*g*lchoose(TotalSize,kSize),
                    BICq=log(n)*kSize - 2*kSize*log(q/(1-q)),
                    CV = 0)
        Criterion <- Neg2LL + penalty
        AllModels[i+1] <- Criterion
        logLAllModels[i+1] <- L0
        if (Criterion<BestModelIC){
                BestModel <-ans
                BestModelIC <- Criterion
                ICLabel3 <- paste("Sample Size: ",n,", Inputs: ",kSize, ", Penalty:",penalty, sep="")
                }
        } #end: for (i in 1:M)
    OMODELS <- order(AllModels)
    BestModels <- matrix(logical(p*TopModels2), ncol=p)
    if (!SFQ)
        dimnames(BestModels)[[2]]<-as.list(names(Xy)[-ncol(Xy),drop=FALSE])
    else
        dimnames(BestModels)[[2]]<-as.list(names(X))
    BestModels<-as.data.frame(BestModels)
    if (M > 1 && TopModels2 > 1) {
         for (i in 1:TopModels2){
            vars<-as.logical(rev(to.binary(OMODELS[i]-1, p)))
            BestModels[i,]<-vars
            }
        BestModels<-cbind(as.data.frame(BestModels), Criterion=AllModels[OMODELS[1:TopModels2]])
    }
    else
        BestModels<-NA
#
#Best subsets of size k
    Subsets <- matrix(logical(p*(p+1)), ncol=p)#initially set to FALSE
    FoundQ <- logical(p+1)
    Criterion <- rep(Inf, p+1)
    if (!SFQ)
        dimnames(Subsets)[[2]]<-as.list(names(Xy)[-ncol(Xy),drop=FALSE])
    else
        dimnames(Subsets)[[2]]<-as.list(names(X))
    BestModels<-as.data.frame(BestModels)
    i <- 0
    while (any(!FoundQ)){
        i <- i+1
        vars<-as.logical(rev(to.binary(OMODELS[i]-1, p)))
        nvars <- sum(vars)
        if (nvars > nvmax)
            break
        if (!FoundQ[nvars+1]) {
            FoundQ[nvars+1] <- TRUE
            Subsets[nvars+1,]<-vars
            logL[nvars+1] <- logLAllModels[OMODELS[i]]
            Criterion[nvars+1] <- AllModels[OMODELS[i]]
            }
        }
    ICLabel2 <- paste(IC, "-", BestModelIC, sep="")
    Subsets <- cbind(Intercept=rep(TRUE,p+1), Subsets) #as from 'leaps'
    indBest <- order(Criterion)[1]
    bestset <- Subsets[indBest,]
}
#OR do we use leaps-and-bounds algorithm
else { 
# 'regsubsets' in leaps package
    if (p > 1)  {#leaps-and-bounds algorithm used in regsubsets when method="exhaustive"
        really.big<-ifelse(P>50, TRUE, FALSE)
#if X contains any factor variables (only 2 levels allowed), we need to convert to 0-1 variables
    if (any(FactorsQ)){
        BinaryCategoricalQ <- TRUE
        indFactors <- (1:p)[FactorsQ]
        for (i in 1:length(indFactors)) 
            X[,indFactors[i]] <- c(unclass(X[,indFactors[i]])-1)
    }
if (BinaryCategoricalQ)
    cat("Note: binary categorical variables converted to 0-1 so 'leaps' could be used.",fill=TRUE)
#
        if (is.null(weights))
            WEIGHTS<-rep(1,n)
        out <- summary(regsubsets(x = X, y = y, intercept=IncludeInterceptQ, nbest=TopModels2, nvmax=NVMAX,
            method=method, really.big=really.big, weights=WEIGHTS, force.in=NULL))
#Subsets has p rows indicating the presence/absence of each variable
#for best subset of size k=1,...,p. If intercept=TRUE, there are p+1 columns.
#The null model, no variables selected is not included yet.
        Subsets <- out$which
        dimnames(Subsets)[[1]]<-1:nrow(Subsets)
        RSS <- out$rss
#Now we include the case k=0, no variables selected.
        if (IncludeInterceptQ)
            RSS <- c(sum((y-mean(y))^2), RSS)
        else {
            RSS <- c(sum(y^2), RSS)
            #Note: 'null model' included in all models like intercept
            Subsets<-cbind(rep(TRUE,nrow(Subsets)), Subsets)
            dimnames(Subsets)[[2]][1]<-"Null"
            }
        logL <- (-n/2)*log(RSS/n) #logL used for qk
        Subsets <- rbind(c(TRUE,rep(FALSE,p)), Subsets)
        dimnames(Subsets)[[1]]<-0:(nrow(Subsets)-1)
        if (TopModels2 > 1) {
            logL <- (-n/2)*log(RSS/n)
            Neg2LL <- -2*logL
            kSize <- rowSums(Subsets)-1 
            penalty <- switch(IC,
                    AIC=2*kSize,
                    BIC=log(n)*kSize,
                    BICg=log(n)*kSize + 2*g*lchoose(p,kSize),
                    BICq=log(n)*kSize - 2*kSize*log(q/(1-q)))
            Criterion <- Neg2LL + penalty
            indTop <- order(Criterion)[1:TopModels2]
            BestModels <- Subsets[indTop,]
            BestModels <- BestModels[,-1]
            BestModels<-cbind(as.data.frame(BestModels), Criterion=Criterion[indTop])
            row.names(BestModels)<-1:TopModels2
            #indBestSize - the best model of each size, k=0,1,...,p
            indBestSize <- match(0:p, kSize)
            Subsets <- Subsets[indBestSize,]
            dimnames(Subsets)[[1]]<-0:p
            RSS <- RSS[indBestSize]
            logL <- logL[indBestSize]
            }
    }
    else { #p=1. gracefully handle this exception.
        Subsets <- matrix(c(TRUE,TRUE), nrow=1, ncol=2)
        if (IncludeInterceptQ)
            colnames(Subsets)<-c("Intercept", colnames(X))
        else
            colnames(Subsets)<-c("Null", colnames(X))
        RSS <- sum((lsfit(X,y)$residuals)^2)
        }
    }
#
#At this point, Subsets, has been obtained either via 'leaps' or enumeration
kSize <- rowSums(Subsets)-1
#Cross-validation methods
if (CVQ){
    CVout <- matrix(numeric(2*P),ncol=2)
    SubsetsAdj<-Subsets[,-1,drop=FALSE]
    for (j in 1:P) {
        if (IC=="CV") {
            CVout[j,] <- switch(CVMethod,
                HTF =  CVHTF(X[,SubsetsAdj[j,],drop=FALSE], y, K=K, REP=REP, family=family, ...),
                DH = CVDH( X[,SubsetsAdj[j,],drop=FALSE], y, K=K, REP=REP),
                d = CVd(X[,SubsetsAdj[j,],drop=FALSE], y, d=K, REP=REP, family=family, ...))
            }
        else #LOOCV only other possibity here
            CVout[j,] <- LOOCV(X[,SubsetsAdj[j,],drop=FALSE], y)
        }
    indBest <- ifelse(CVMethod=="HTF", oneSdRule(CVout),which.min(CVout[,1]))
    }
#Information criterion methods. Only for leaps do we need these.
else if (LEAPSQ) {
#Special case: q=0 or q=1
    if (IC=="BICq" && (q==0||q==1)){
        if (q == 0) {
            indBest<-1
            Criterion<-log(0)
            }
        else {
            indBest<-nrow(Subsets)
            Criterion<-log(1/0)
            }
        }
    else {
#obtain the particular information criterion
        Neg2LL <- -2*logL 
        penalty <- switch(IC,
                    AIC=2*kSize,
                    BIC=log(n)*kSize,
                    BICg=log(n)*kSize + 2*g*lchoose(p,kSize),
                    BICq=log(n)*kSize - 2*kSize*log(q/(1-q)))
        Criterion <- Neg2LL + penalty
        indBest <- order(Criterion)[1]
    }
} #end of 'if (LEAPSQ)...'
#Fit 'bestmodel'
if (!glmQ || CVQ){
    bestset <- Subsets[indBest,]
    if (indBest == 1)
        if (gaussianQ) {
            if (IncludeInterceptQ)
                BestModel <- lm(y~1, ...)
            else
                BestModel <- lm(y~-1, ...)
            }
        else {
            if (IncludeInterceptQ)
                BestModel <- glm(y~1, family=family, ...)
            else
                BestModel <- glm(y~-1, family=family, ...)
            }
    else
        if (gaussianQ){
            if (IncludeInterceptQ)
                BestModel <-lm(y~., data=data.frame(Xy[,c(bestset[-1],FALSE),drop=FALSE],y=y), ...)
            else
                BestModel <-lm(y~-1+., data=data.frame(Xy[,c(bestset[-1],FALSE),drop=FALSE],y=y), ...)
            }
        else {
            if (IncludeInterceptQ)
                BestModel <-glm(y~., family=family, data=data.frame(Xy[,c(bestset[-1],FALSE),drop=FALSE],y=y), ...)
            else
                BestModel <-glm(y~-1+., family=family, data=data.frame(Xy[,c(bestset[-1],FALSE),drop=FALSE],y=y), ...)
            }
    }
#
ModelReport$Bestk <- sum(as.numeric(bestset))-1
Subsets <- switch(IC,
    AIC  =cbind(as.data.frame(Subsets), logLikelihood=logL, AIC=Criterion),
    BIC  =cbind(as.data.frame(Subsets), logLikelihood=logL, BIC=Criterion),
    BICg =cbind(as.data.frame(Subsets), logLikelihood=logL, BICg=Criterion),
    BICq =cbind(as.data.frame(Subsets), logLikelihood=logL, BICq=Criterion),
    LOOCV =cbind(as.data.frame(Subsets), logLikelihood=logL, LOOCV=CVout[,1]),
    CV   =switch(CVMethod,
            HTF=cbind(as.data.frame(Subsets), logLikelihood=logL, CV=CVout[,1], sdCV=CVout[,2]),
            DH =cbind(as.data.frame(Subsets), logLikelihood=logL, CV=CVout[,1]),
            d =cbind(as.data.frame(Subsets),  logLikelihood=logL, CV=CVout[,1])))
# 
#label rows to indicate best model
row.names(Subsets)<-paste(0:p, ifelse(indBest==1:P,"*"," "),sep="")
#Obtain q tuning table
#q tuning table not available with some options:
if ((CategoricalQ && !BinaryCategoricalQ) || p==1 || NVMAX != p )
    Bestq<-qk<-NA #maybe someday this will be improved
else {
    q12 <- sapply(2:(P-1),function(k){
        i1<-1:(k-1) 
        i2<-(k+1):P
        0.5*log(n)-c(min((logL[i1]-logL[k])/(i1-k)), max((logL[i2]-logL[k])/(i2-k)))
        }   
        )
    q12 <- matrix(q12, ncol=2, byrow=TRUE)
    q12 <- rbind(    
        c(-Inf, 
        0.5*log(n)-max((logL[2:P]-logL[1])/(2:P-1))),
        q12, 
        c(0.5*log(n)-min((logL[1:(P-1)]-logL[P])/(1:(P-1)-P)),
        Inf)
        )
    indNum<-!apply(is.na(q12), MARGIN=1, FUN=any)
    q12<-q12[indNum,]
#numerically better, Changjiang email June 12
    q12<- 1-1/(1+exp(q12))
    #q12 <- exp(q12)/(1+exp(q12))
    #correction for upper limit
    q12[nrow(q12),2] <- 1  
    ind <- q12[,1]<q12[,2]
    QK <- cbind(logL, q12, (0:p))
    qk <- QK[ind,]
#
#find optimal q
    level<-qLevel
    ULBound<- qchisq(level,1)
    ks<- 0:p
    D12<- sapply(2:(P-1),function(k){
        i1<- 1:(k-1)
        i2<- (k+1):P
        2*c(max((logL[i2]-logL[k])/(i2-k)), min((logL[i1]-logL[k])/(i1-k)))
        } )
    D12<- matrix(D12,ncol=2,byrow=TRUE)
    D12<- rbind(2*c(max((logL[2:P]-logL[1])/(2:P-1)),Inf),D12, 2*c(0, min((logL[1:(P-1)]-logL[P])/(1:(P-1)-P))))
    kS1<- max(which(D12[,1]<D12[,2]&D12[,2]>=ULBound))
    DQ<- D12[,2]-D12[,1]
    kS2<- which(DQ==max(DQ[D12[,1]<ULBound]))
#
    Bestq<-matrix(cbind(q12[c(kS1,kS2),],ks[c(kS1,kS2)]), nrow=2)
    dimnames(Bestq) <- list(c("BICq1","BICq2"), c("q1","q2","selected k"))
    indq <- match(kSize[indBest], qk[,4])
    if (is.na(indq))
        ICLabel <- paste(ICLabel, "\nNo BICq equivalent", sep="")
    else {
        BICqLab<-paste("BICq equivalent for q in (", qk[indq,2],", ",qk[indq,3],")",sep="")
        ICLabel <- paste(ICLabel, "\n", BICqLab, sep="")
    }
    colnames(qk) <-c("LogL", "q1", "q2", "k")
    rownames(qk) <-NULL
}
Subsets<-Subsets[is.finite(Subsets[,P+1]),]
if (!IncludeInterceptQ)
    Subsets<-Subsets[,-1]
out <- list(BestModel=BestModel, BestModels=BestModels, Bestq=Bestq, qTable=qk, Subsets=Subsets,    
            Title=ICLabel, ModelReport=ModelReport)
class(out) <- "bestglm"
out
}
