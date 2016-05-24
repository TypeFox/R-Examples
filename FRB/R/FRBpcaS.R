FRBpcaS <- function(Y,...) UseMethod("FRBpcaS")


FRBpcaS.formula <- function (formula, data = NULL, ...) 
{
.check_vars_numeric <- function (mf) 
{
    mt <- attr(mf, "terms")
    mterms <- attr(mt, "factors")
    mterms <- rownames(mterms)[apply(mterms, 1, any)]
    any(sapply(mterms, function(x) is.factor(mf[, x]) || !is.numeric(mf[, 
        x])))
}
#--------------------------------------------------------------------------
# Main function

    cl <- match.call()
    mt <- terms(formula, data = data)
    if (attr(mt, "response") > 0L) 
        stop("response not allowed in formula")
    mf <- match.call(expand.dots = FALSE)
    mf$... <- NULL
    mf[[1L]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    if (.check_vars_numeric(mf)) 
        stop("PCA applies only to numerical variables")
    na.act <- attr(mf, "na.action")
    mt <- attr(mf, "terms")
    attr(mt, "intercept") <- 0L
    x <- model.matrix(mt, mf)
    res <- FRBpcaS.default(x, ...)
    cl <- match.call()
    cl[[1]] <- as.name("FRBpcaS")
    res$call <- cl
    if (!is.null(na.act)) res$na.action <- na.act
    return(res)
}

FRBpcaS.default <- function(Y, R=999, bdp=0.5, conf=0.95, control=Scontrol(...), na.action=na.omit, ...)
{
# performs PCA based on the multivariate S estimate of shape, with
# fast and robust bootstrap
#
# calls: Sest_loccov(), Sboot_loccov(), Seinfs_pca()
#
# INPUT :
#   Y : (n x q) data
#   R : number of bootstrap samples
#   bdp : breakdown point of S-estimate (determines tuning parameters)
#   conf : confidence level for bootstrap intervals
# OUTPUT :
#   res$est : (list) result of Sest_loccov()
#   res$bootest : (list) result of Sboot_loccov()
#   res$shape : (q x q) S-estimate of the shape matrix
#   res$eigval : (q x 1) eigenvalues of S shape
#   res$eigvec : (q*q x 1) eigenvectors of S-shape
#   res$pvar : (q-1 x 1) percentages of variance for S eigenvalues
#   res$eigval.boot : (q x R) eigenvalues of S shape
#   res$eigvec.boot : (q*q x R) eigenvectors of S-shape
#   res$pvar.boot : (q-1 x R) percentages of variance for S eigenvalues
#   res$eigval.SE : (q x 1) bootstrap variance for S eigenvalues
#   res$eigvec.SE : (q x q) bootstrap variance for S eigenvectors
#   res$pvar.SE : (q-1 x 1) bootstrap variance for percentage of
#                                               variance for S-eigenvalues
#   res$avgangle : (q x 1) average angles between bootstrap eigenvectors
#                                           and original S eigenvectors
#   res$eigval.CI.bca : (q x 2) % BCa intervals for S eigenvalues
#   res$eigvec.CI.bca : (q*q x 2) % BCa intervals for S eigenvectors
#   res$pvar.CI.bca : (q-1 x 2) % BCa intervals for percentage of
#                                               variance for S-eigenvalues
#   res$pvar.CIone.bca : (q-1 x 1) % one-sided BCa intervals for 
#               percentage of variance for S-eigenvalues ([-infty upper]) 
#   res$eigval.CI.basic : (q x 2) % basic bootstrap intervals for S eigenvalues
#   res$eigvec.CI.basic : (q*q x 2) % basic bootstrap intervals for S eigenvectors
#   res$pvar.CI.basic : (q-1 x 2) % basic bootstrap intervals for percentage of
#                                               variance for S-eigenvalues
#   res$pvar.CIone.basic : (q-1 x 1) % one-sided basic bootstrap intervals for 
#               percentage of variance for S-eigenvalues ([-infty upper]) 
#   res$failedsamples : number of discarded bootstrap samples due to
#                                   non-positive definiteness of shape
#                   (results from 'discarded' bootstrap samples are omitted
#            such that actual number of recalculations can be lower than R) 
#   REMARK: whenever matrices contain eigenvalue-related values, the ordering is "DECREASING"
#                       (e.g. first column of eigenvectors contains vector corresponding to largest eigenvalue)                
#                       (or first row in case of the matrix of bootstrapped values)


# --------------------------------------------------------------------

vecop <- function(mat) {
# performs vec-operation (stacks colums of a matrix into column-vector)

nr <- nrow(mat)
nc <- ncol(mat)

vecmat <- rep(0,nr*nc)
for (col in 1:nc) {
    startindex <- (col-1)*nr+1
    vecmat[startindex:(startindex+nr-1)] <- mat[,col]
}
return(vecmat)
}

# --------------------------------------------------------------------

reconvec <- function(vec,ncol) {
# reconstructs vecop'd matrix

lcol <- length(vec)/ncol
rec <- matrix(0,lcol,ncol)
for (i in 1:ncol)
    rec[,i] <- vec[((i-1)*lcol+1):(i*lcol)]

return(rec)
}

# --------------------------------------------------------------------
# -                        main function                             -
# --------------------------------------------------------------------

Y <- as.matrix(Y)
Y=na.action(Y)

n <- nrow(Y)
q <- ncol(Y)
if (q <= 1L) stop("at least two variables needed for PCA")
if (n < q) stop("For PCA the number of observations cannot be smaller than the 
number of variables")

#bdp <- .5
dimens <- q + q*q

# compute S-estimates of location and shape
Sests <- Sest_loccov(Y, bdp=bdp, control=control)

SGamma <- Sests$Gamma
SSigma <- SGamma * Sests$scale^2

# compute corresponding original S-eigenvalue/eigenvector estimates
eigresGammaS <- eigen(SGamma)
IXGamma <- order(eigresGammaS$values, decreasing=TRUE)
SeigvalsGamma <- eigresGammaS$values[IXGamma]
Seigvecs <- eigresGammaS$vectors[,IXGamma]

majorloading <- rep(0,q)
for (k in 1:q) {
    IX <- order(abs(Seigvecs[,k]))
    majorloading[k] <- IX[q]    
    Seigvecs[,k] <- sign(Seigvecs[IX[q],k])*Seigvecs[,k]
}  
vecSeigvecs <- vecop(Seigvecs)

Svarperc <- rep(0,q-1)
for (k in 1:(q-1)) {
    Svarperc[k] <- sum(SeigvalsGamma[1:k])/sum(SeigvalsGamma)
}

# compute bootstrapped S-estimates of location and shape
bootres <- Sboot_loccov(Y, R, ests=Sests)

#################################################################################
# now take the bootstrapped shape estimates and compute their eigenvalues/vectors

booteigvalsGamma <- matrix(0,q,R)
booteigvecs <- matrix(0,q*q,R)
bootpercvars <- matrix(0,q-1,R)
bootangles <- matrix(0,q,R)
bootsampleOK <- rep(1,R)
bootsampleOKreally <- rep(1,R)

for (r in 1:R) {
    correctedSSigmast <- reconvec(bootres$centered[(q+1):dimens,r],q) + SSigma
    detSigma <- det(correctedSSigmast)
    if (detSigma<0) {
        bootsampleOK[r] <- 0
        correctedSSigmast <- make.positive.definite(correctedSSigmast)
        detSigma <- prod(eigen(correctedSSigmast)$values) # to be safe
        if (is.complex(detSigma) | Re(detSigma)<=0) {bootsampleOKreally[r] <- 0; next}
    }
    correctedSGammast <- detSigma^(-1/q) * correctedSSigmast
    eigresGammast <- eigen(correctedSGammast)
IXGammast <- order(eigresGammast$values, decreasing=TRUE)    
    eigenvaluesGammast <- eigresGammast$values[IXGammast]
    if (any(is.complex(eigenvaluesGammast))){bootsampleOKreally[r] <- 0; next}
    if (any(eigenvaluesGammast<0))	  {bootsampleOKreally[r] <- 0; next} 
#        eigenvaluesGammast <- pmax(0, eigenvaluesGammast)
#        bootsampleOK[r] <- 0
#    }
    eigenvectorsGammast <- eigresGammast$vectors[,IXGammast]
    for (k in 1:q) {
        eigenvectorsGammast[,k] <- sign(eigenvectorsGammast[majorloading[k],k]) * eigenvectorsGammast[,k]
    }    
 
    percvarGammast <- rep(0,q-1)
    for (k in 1:(q-1)) {
        percvarGammast[k] <- sum(eigenvaluesGammast[1:k])/sum(eigenvaluesGammast)
    }

    Svecs <- eigenvectorsGammast
    for (k in 1:q) {
        bootangles[k,r] <- acos(min(abs(t(Svecs[,k]) %*% Seigvecs[,k]/sqrt(t(Svecs[,k]) %*% Svecs[,k])/sqrt(t(Seigvecs[,k])%*%Seigvecs[,k])),1))
    }
    booteigvalsGamma[,r] <- eigenvaluesGammast
    booteigvecs[,r] <- vecop(eigenvectorsGammast)
    bootpercvars[,r] <- percvarGammast
}

# we discard bootstrap samples with non-positive definite covariance S-estimate
bootindicesOK <- (1:R)[bootsampleOK==1]
nfailed <- R - length(bootindicesOK)
# ... except if there were too many of those...
if (nfailed > 0.75*R) {
    warning("more than 75% of bootstrapped shape matrices was non-positive definite; 
    they were used anyway in case make.positive.definite was succesful")
    bootindicesOK <- (1:R)[bootsampleOKreally==1]
    nfailed <- R - length(bootindicesOK)
}

bootangles <- bootangles[,bootindicesOK]  
booteigvalsGamma <- booteigvalsGamma[,bootindicesOK]
bootpercvars <- bootpercvars[,bootindicesOK,drop=FALSE]
booteigvecs <- booteigvecs[,bootindicesOK]
                                    
avgangle <- apply(bootangles,1,mean)

# compute bootstrap estimates of variance
eigGvariances <- apply(booteigvalsGamma,1,var)
pvarvariances <- apply(bootpercvars,1,var)
eigvecvariances <- apply(booteigvecs,1,var)

# sort bootstrap recalculations for constructing intervals
sortedeigG <- t(apply(booteigvalsGamma,1,sort))
sortedpvar <- t(apply(bootpercvars,1,sort))
sortedeigvec <- t(apply(booteigvecs,1,sort))

##################################################################################
# compute BCa confidence limits

# empirical inlfuences for computing a in BCa intervals, based on IF(S)
Einf <- Seinfs_pca(Y, ests=Sests)
eigGinflE <- Einf$eigs
eigvecinflE <- Einf$eigvec
pvarinflE <- Einf$varperc

# set Rok equal to the actual number of OK bootstrap samples
Rok <- length(bootindicesOK)

normquan <- qnorm(1 - (1 - conf)/2)
normquanone <- qnorm(conf)
basicindexlow <- floor((1 - (1 - conf)/2) * Rok)
basicindexhigh <- ceiling((1 - conf)/2 * Rok)
basicindexhighone <- ceiling((1 - conf) * Rok)

eigGCIbca <- matrix(0,q,2)
eigGCIbasic <- matrix(0,q,2)
for (i in 1:q) {
    nofless <- length(sortedeigG[i,sortedeigG[i,] <= SeigvalsGamma[i]])
    w <- qnorm((nofless+1)/(Rok+2))
    a <- 1/6 * sum(eigGinflE[,i]^3) / (sum(eigGinflE[,i]^2)^(3/2))
    alphatildelow <- pnorm(w+(w-normquan)/(1-a*(w-normquan)))
    alphatildehigh <- pnorm(w+(w+normquan)/(1-a*(w+normquan)))
    indexlow <- min(max((Rok+1)*alphatildelow,1), Rok)
    indexhigh <- max(min((Rok+1)*alphatildehigh,Rok),1)
    eigGCIbca[i,1] <- sortedeigG[i,round(indexlow)]
    eigGCIbca[i,2] <- sortedeigG[i,round(indexhigh)]
}
eigGCIbasic[,1] <- 2*SeigvalsGamma - sortedeigG[,basicindexlow]
eigGCIbasic[,2] <- 2*SeigvalsGamma - sortedeigG[,basicindexhigh]

eigvecCIbca <- matrix(0,q*q,2)
eigvecCIbasic <- matrix(0,q*q,2)
for (i in 1:(q*q)) {
    nofless <- length(sortedeigvec[i,sortedeigvec[i,] <= vecSeigvecs[i]])
    w <- qnorm((nofless+1)/(Rok+2))
    a <- 1/6 * sum(eigvecinflE[,i]^3) / (sum(eigvecinflE[,i]^2)^(3/2))
    alphatildelow <- pnorm(w+(w-normquan)/(1-a*(w-normquan)))
    alphatildehigh <- pnorm(w+(w+normquan)/(1-a*(w+normquan)))
    indexlow <- min(max((Rok+1)*alphatildelow,1), Rok)
    indexhigh <- max(min((Rok+1)*alphatildehigh,Rok),1)
    eigvecCIbca[i,1] <- sortedeigvec[i,round(indexlow)]
    eigvecCIbca[i,2] <- sortedeigvec[i,round(indexhigh)]
}
eigvecCIbasic[,1] <- 2*vecSeigvecs - sortedeigvec[,basicindexlow]
eigvecCIbasic[,2] <- 2*vecSeigvecs - sortedeigvec[,basicindexhigh]

pvarCIbca <- matrix(0,q-1,2)
pvarCIbcaone <- matrix(0,q-1,1)
pvarCIbasic <- matrix(0,q-1,2)
pvarCIbasicone <- matrix(0,q-1,1)
for (i in 1:(q-1)) {
    nofless <- length(sortedpvar[i,sortedpvar[i,] <= Svarperc[i]])
    w <- qnorm((nofless+1)/(Rok+2))
    a <- 1/6 * sum(pvarinflE[,i]^3) / (sum(pvarinflE[,i]^2)^(3/2))
    alphatildelow <- pnorm(w+(w-normquan)/(1-a*(w-normquan)))
    alphatildehigh <- pnorm(w+(w+normquan)/(1-a*(w+normquan)))
    indexlow <- min(max((Rok+1)*alphatildelow,1), Rok)
    indexhigh <- max(min((Rok+1)*alphatildehigh,Rok),1)
    pvarCIbca[i,1] <- sortedpvar[i,round(indexlow)]
    pvarCIbca[i,2] <- sortedpvar[i,round(indexhigh)]

    # one-sided interval
    alphatildehigh <- pnorm(w+(w+normquanone)/(1-a*(w+normquanone)))
    indexhigh <- max(min((Rok+1)*alphatildehigh,Rok),1)
    pvarCIbcaone[i] <- sortedpvar[i,round(indexhigh)]
}
pvarCIbasic[,1] <- 2*Svarperc - sortedpvar[,basicindexlow]
pvarCIbasic[,2] <- 2*Svarperc - sortedpvar[,basicindexhigh]
pvarCIbasicone <- 2*Svarperc - sortedpvar[,basicindexhighone]

####################################################################################

method <- paste("PCA based on multivariate S-estimates (breakdown point = ", bdp, ")", sep="")

z <- list(est=Sests, bootest=bootres, shape=SGamma, eigval=SeigvalsGamma, eigvec=Seigvecs,
        pvar=Svarperc, eigval.boot=booteigvalsGamma, eigvec.boot=booteigvecs, pvar.boot=bootpercvars, 
        eigval.SE=sqrt(eigGvariances), eigvec.SE=sqrt(reconvec(eigvecvariances,q)), pvar.SE=sqrt(pvarvariances), 
        angles=bootangles, avgangle=avgangle, eigval.CI.bca=eigGCIbca, eigvec.CI.bca=eigvecCIbca, pvar.CI.bca=pvarCIbca, pvar.CIone.bca=pvarCIbcaone, 
        eigval.CI.basic=eigGCIbasic, eigvec.CI.basic=eigvecCIbasic, pvar.CI.basic=pvarCIbasic, pvar.CIone.basic=pvarCIbasicone, 
        failedsamples=nfailed, conf=conf, method=method, w=Sests$w, outFlag=Sests$outFlag, Y=Y)
        
class(z) <- "FRBpca"    

return(z)    
        
         
}
