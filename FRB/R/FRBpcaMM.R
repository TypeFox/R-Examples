FRBpcaMM <- function(Y,...) UseMethod("FRBpcaMM")


FRBpcaMM.formula <- function (formula, data = NULL, ...) 
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
    res <- FRBpcaMM.default(x, ...)
  	cl <- match.call()
    cl[[1]] <- as.name("FRBpcaMM")
    res$call <- cl
    if (!is.null(na.act)) res$na.action <- na.act
    return(res)
}



FRBpcaMM.default <- function(Y, R=999, conf=0.95, control=MMcontrol(...),na.action=na.omit, ...)
{
# performs PCA based on the multivariate MM estimate of shape, with
# fast and robust bootstrap
#
# calls: MMest_loccov(), MMboot_loccov(), MMeinfs_pca()
#
# INPUT :
#   Y : (n x q) data
#   R : number of bootstrap samples
#   conf : confidence level for bootstrap intervals
# OUTPUT :
#   res$est : (list) result of MMest_loccov()
#   res$bootest : (list) result of MMboot_loccov()
#   res$shape : (q x q) MM-estimate of the shape matrix
#   res$eigval : (q x 1) eigenvalues of MM shape
#   res$eigvec : (q x q) eigenvectors of MM-shape
#   res$pvar : (q-1 x 1) percentages of variance for MM eigenvalues
#   res$eigval.boot : (q x R) eigenvalues of MM shape
#   res$eigvec.boot : (q*q x R) eigenvectors of MM-shape
#   res$pvar.boot : (q-1 x R) percentages of variance for MM eigenvalues
#   res$eigval.SE : (q x 1) bootstrap standard error for MM eigenvalues
#   res$eigvec.SE : (q x q) bootstrap standard error for MM eigenvectors
#   res$pvar.SE : (q-1 x 1) bootstrap standard error for percentage of
#                                               variance for MM-eigenvalues
#   res$avgangle : (q x 1) average angles between bootstrap eigenvectors
#                                           and original MM eigenvectors
#   res$eigval.CI.bca : (q x 2) % BCa intervals for MM eigenvalues
#   res$eigvec.CI.bca : (q*q x 2) % BCa intervals for MM eigenvectors
#   res$pvar.CI.bca : (q-1 x 2) % BCa intervals for percentage of
#                                               variance for MM-eigenvalues
#   res$pvar.CIone.bca : (q-1 x 1) % one-sided BCa intervals for 
#               percentage of variance for MM-eigenvalues ([-infty upper]) 
#   res$eigval.CI.basic : (q x 2) % basic bootstrap intervals for MM eigenvalues
#   res$eigvec.CI.basic : (q*q x 2) % basic bootstrap intervals for MM eigenvectors
#   res$pvar.CI.basic : (q-1 x 2) % basic bootstrap intervals for percentage of
#                                               variance for MM-eigenvalues
#   res$pvar.CIone.basic : (q-1 x 1) % one-sided basic bootstrap intervals for 
#               percentage of variance for MM-eigenvalues ([-infty upper]) 
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

dimens <- q + q*q

# compute MM-estimates of location and shape
control$shapeEff <- TRUE # suppressing the input!
MMests <- MMest_loccov(Y, control=control)
MMGamma <- MMests$Gamma

# compute corresponding original MM- eigenvalue/eigenvector estimates
eigresGammaMM <- eigen(MMGamma)
IXGamma <- order(eigresGammaMM$values, decreasing=TRUE)
MMeigvalsGamma <- eigresGammaMM$values[IXGamma]
MMeigvecs <- eigresGammaMM$vectors[,IXGamma]
                                                          
majorloading <- rep(0,q)
for (k in 1:q) {
    IX <- order(abs(MMeigvecs[,k]))
    majorloading[k] <- IX[q]    
    MMeigvecs[,k] <- sign(MMeigvecs[IX[q],k])*MMeigvecs[,k]
}  
vecMMeigvecs <- vecop(MMeigvecs)

MMvarperc <- rep(0,q-1)
for (k in 1:(q-1)) {
    MMvarperc[k] <- sum(MMeigvalsGamma[1:k])/sum(MMeigvalsGamma)
}

# compute bootstrapped MM-estimates of location and shape
bootres <- MMboot_loccov(Y, R, ests=MMests)

##################################################################################
# now take the bootstrapped shape estimates and compute their eigenvalues/vectors

booteigvalsGamma <- matrix(0,q,R)
booteigvecs <- matrix(0,q*q,R)
bootpercvars <- matrix(0,q-1,R)
bootangles <- matrix(0,q,R)
bootsampleOK <- rep(1,R)

for (r in 1:R) {
    correctedMMGammast <- reconvec(bootres$centered[(q+1):dimens,r],q) + MMGamma
    eigresGammast <- eigen(correctedMMGammast)
    IXGammast <- order(eigresGammast$values, decreasing=TRUE)    
    eigenvaluesGammast <- eigresGammast$values[IXGammast]
    if (any(eigenvaluesGammast<0)) { 
        eigenvaluesGammast <- pmax(0, eigenvaluesGammast)
        bootsampleOK[r] <- 0
    }
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
        bootangles[k,r] <- acos(min(abs(t(Svecs[,k]) %*% MMeigvecs[,k]/sqrt(t(Svecs[,k]) %*% Svecs[,k])/sqrt(t(MMeigvecs[,k])%*%MMeigvecs[,k])),1))
    }
    booteigvalsGamma[,r] <- eigenvaluesGammast
    booteigvecs[,r] <- vecop(eigenvectorsGammast)
    bootpercvars[,r] <- percvarGammast
}


# we discard bootstrap samples with non-positive definite covariance MM-estimate
bootindicesOK <- (1:R)[bootsampleOK==1]
nfailed <- R - length(bootindicesOK)
if (nfailed > 0.75*R) {
    warning("more than 75% of bootstrapped shape matrices was non-positive definite; they were used anyway after eigenvalue adjustment")
    bootindicesOK <- 1:R
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

# empirical inlfuences for computing a in BCa intervals, based on IF(MM)
EinfMM <- MMeinfs_pca(Y, ests=MMests)
eigGinflE <- EinfMM$eigs
eigvecinflE <- EinfMM$eigvec
pvarinflE <- EinfMM$varperc

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
    nofless <- length(sortedeigG[i,sortedeigG[i,]<=MMeigvalsGamma[i]])
    w <- qnorm((nofless+1)/(Rok+2))
    a <- 1/6 * sum(eigGinflE[,i]^3) / (sum(eigGinflE[,i]^2)^(3/2))
    alphatildelow <- pnorm(w+(w-normquan)/(1-a*(w-normquan)))
    alphatildehigh <- pnorm(w+(w+normquan)/(1-a*(w+normquan)))
    indexlow <- min(max((Rok+1)*alphatildelow,1), Rok)
    indexhigh <- max(min((Rok+1)*alphatildehigh,Rok),1)
    eigGCIbca[i,1] <- sortedeigG[i,round(indexlow)]
    eigGCIbca[i,2] <- sortedeigG[i,round(indexhigh)]
}
eigGCIbasic[,1] <- 2*MMeigvalsGamma - sortedeigG[,basicindexlow]
eigGCIbasic[,2] <- 2*MMeigvalsGamma - sortedeigG[,basicindexhigh]

eigvecCIbca <- matrix(0,q*q,2)
eigvecCIbasic <- matrix(0,q*q,2)
for (i in 1:(q*q)) {
    nofless <- length(sortedeigvec[i,sortedeigvec[i,]<=vecMMeigvecs[i]])
    w <- qnorm((nofless+1)/(Rok+2))
    a <- 1/6 * sum(eigvecinflE[,i]^3) / (sum(eigvecinflE[,i]^2)^(3/2))
    alphatildelow <- pnorm(w+(w-normquan)/(1-a*(w-normquan)))
    alphatildehigh <- pnorm(w+(w+normquan)/(1-a*(w+normquan)))
    indexlow <- min(max((Rok+1)*alphatildelow,1), Rok)
    indexhigh <- max(min((Rok+1)*alphatildehigh,Rok),1)
    eigvecCIbca[i,1] <- sortedeigvec[i,round(indexlow)]
    eigvecCIbca[i,2] <- sortedeigvec[i,round(indexhigh)]
}
eigvecCIbasic[,1] <- 2*vecMMeigvecs - sortedeigvec[,basicindexlow]
eigvecCIbasic[,2] <- 2*vecMMeigvecs - sortedeigvec[,basicindexhigh]

pvarCIbca <- matrix(0,q-1,2)
pvarCIbcaone <- matrix(0,q-1,1)
pvarCIbasic <- matrix(0,q-1,2)
pvarCIbasicone <- matrix(0,q-1,1)
for (i in 1:(q-1)) {
    nofless <- length(sortedpvar[i,sortedpvar[i,]<=MMvarperc[i]])
    w <- qnorm((nofless+1) /(Rok+2))
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
pvarCIbasic[,1] <- 2*MMvarperc - sortedpvar[,basicindexlow]
pvarCIbasic[,2] <- 2*MMvarperc - sortedpvar[,basicindexhigh]
pvarCIbasicone <- 2*MMvarperc - sortedpvar[,basicindexhighone]

####################################################################################

method <- paste("PCA based on multivariate MM-estimates (bdp = ", control$bdp, ", eff = ", control$eff, ")", sep="")

z <- list(est=MMests, bootest=bootres, shape=MMGamma, eigval=MMeigvalsGamma, eigvec=MMeigvecs, pvar=MMvarperc, eigval.boot=booteigvalsGamma, 
        eigvec.boot=booteigvecs, pvar.boot =bootpercvars, eigval.SE=sqrt(eigGvariances), eigvec.SE=sqrt(reconvec(eigvecvariances,q)),
        pvar.SE=sqrt(pvarvariances), angles=bootangles, avgangle=avgangle, eigval.CI.bca=eigGCIbca, eigvec.CI.bca=eigvecCIbca, pvar.CI.bca=pvarCIbca, 
        pvar.CIone.bca=pvarCIbcaone, eigval.CI.basic=eigGCIbasic, eigvec.CI.basic=eigvecCIbasic, pvar.CI.basic=pvarCIbasic, 
        pvar.CIone.basic=pvarCIbasicone, failedsamples=nfailed, conf=conf, method=method, w=MMests$w, outFlag=MMests$outFlag, Y=Y)

class(z) <- "FRBpca"    

return(z)    
}
