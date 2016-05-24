## Rao's score test statistic with the most probable copy numbers

##' Calculates the score test statistics with the most probable CNV.
##'
##' @title Score test with the most probable CNV
##' @param envirX The matrix of environmental variables. The intercept should be included if it's needed.
##' @param fam The FAM file which follows the format defined in PLINK.
##' @param clusRes The clustering group which is signed to each individual.
##' @param alpha The estimated parameters for environmental variables under null hypothesis. This value can be calculated by using function \code{\link{AssoTestProc}}. 
##' @param phi The matrix of correlation between individuals.
##' @param sig2g The estimated standard error for polygenic effect under null hypothesis. This value can be calculated by using function \code{\link{AssoTestProc}}. 
##' @param sig2 The estimated standard error for environmental effect under null hypothesis. This value can be calculated by using function \code{\link{AssoTestProc}}. 
##' @return It returns the statistic value and pvalue of the score test. 
##' \item{STEs}{The statistic value of score test with the most probable CNV.}
##' \item{STEp}{The pvalue of score test with the most probable CNV.}
##' @author Meiling Liu, Sungho Won 
##' @examples
##' # Fit the data under the assumption that there are 3 clusters
##' asso.fit <- AssoTestProc(signal=signal,fam=fam,envirX=envirX,phi=phi,N=3,varSelection='PC.9')
##' cnv_e <- asso.fit$clusRes
##' alpha <- asso.fit$para$alpha
##' sig2g <- asso.fit$para$sig2g
##' sig2 <- asso.fit$para$sig2
##' STE(envirX=envirX,clusRes=cnv_e,fam=fam,alpha=alpha,phi=phi,sig2g=sig2g,sig2=sig2)
##' @export

STE<- function(envirX,clusRes,fam,alpha,phi,sig2g,sig2){

    envirX <- cbind(1,envirX)
    pheno <- as.numeric(fam[,6])
    S  <- length(pheno)
    Invcov <- solve(sig2g*phi+sig2*diag(S))
    e <- pheno-envirX%*%alpha

    nve <- t(clusRes)%*%Invcov%*%e
    nvn <- t(clusRes)%*%Invcov%*%clusRes
    nvx <- t(clusRes)%*%Invcov%*%envirX
    xvx <- t(envirX)%*%Invcov%*%envirX
    Trs <- nve%*%solve(nvn-nvx%*%solve(xvx)%*%t(nvx))%*%t(nve)
    pv <- 1-pchisq(Trs,1)
    return(list(STEs=Trs,STEp=pv))
}

##  Rao's score test statistic with the probe intensity measurements

##' Calculates the score test statistics with the intensity value.
##'
##' @title Score test with the intensity value
##' @param envirX The matrix of environmental variables. The intercept should be included if it's needed.
##' @param fam The FAM file which follows the format defined in PLINK.
##' @param signal The matrix of intensity measurements. The row names must be consistent with the Individual ID in fam file.
##' @param alpha The estimated parameters for environmental variables under null hypothesis. This value can be calculated by using function \code{\link{AssoTestProc}}. 
##' @param phi The matrix of correlation between individuals.
##' @param sig2g The estimated standard error for polygenic effect under null hypothesis. This value can be calculated by using function \code{\link{AssoTestProc}}. 
##' @param sig2 The estimated standard error for environmental effect under null hypothesis. This value can be calculated by using function \code{\link{AssoTestProc}}. 
##' @return It returns the statistic value and pvalue of the score test. 
##' \item{STIMs}{The statistic value of score test with the intensity value under null hypothesis.}
##' \item{STIMp}{The pvalue of score test with the intensity value under null hypothesis.}
##' \item{df}{The degree of freedom of score test with the intenstiy value under null hypothesis.}
##' @author Meiling Liu, Sungho Won 
##' @examples
##' # Fit the data under the assumption that there are 3 clusters
##' asso.fit <- AssoTestProc(signal=signal,fam=fam,envirX=envirX,phi=phi,N=3,varSelection='PC.9')
##' alpha <- asso.fit$para$alpha
##' sig2g <- asso.fit$para$sig2g
##' sig2 <- asso.fit$para$sig2
##' STIM(envirX=envirX,signal=signal,fam=fam,alpha=alpha,phi=phi,sig2g=sig2g,sig2=sig2)
##' @export

STIM<- function(envirX,signal,fam,alpha,phi,sig2g,sig2){

    

    envirX <- cbind(1,envirX)
    pheno <- as.numeric(fam[,6])
    S  <- length(pheno)
    InvW <- solve(phi)
    Invcov <- solve(sig2g*phi+sig2*diag(S))    
    temp <- scale(signal)
    R <- cor(t(temp))

    N1 <- matrix(1,S,1)
    vz <- Invcov%*%envirX

    s1 <- envirX%*%solve(t(envirX)%*%vz)%*%t(vz)
    InvWN1 <- InvW%*%N1
    s2 <- N1%*%solve(t(N1)%*%InvWN1)%*%t(InvWN1)
    t3 <- t(pheno)%*%t(diag(S)-s1)%*%Invcov%*%(diag(S)-s2)
    t1 <- (diag(S)-s2)%*%R
    t2 <- sum(diag(t1%*%t(diag(S)-s2)%*%Invcov%*%(diag(S)-s1)))
    u <- t3%*%signal
    psi <- t(signal)%*%t(diag(S)-s2)%*%(diag(S)-s2)%*%signal
    v <- t2/sum(diag(t1))*psi
    Trs <- u%*%solve(v)%*%t(u)
    df <- qr(v)$rank
    pv<- 1-pchisq(Trs,df)

    return(list(STIMs=Trs,STIMp=pv,df=df))
}



