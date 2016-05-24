GSboot_multireg <- function(X,Y,R=999,conf=0.95,ests=GSest_multireg(X,Y))
{
# robust bootstrap for multivariate GS-regression (only for the slope)+
# confidence intervals
# INPUT:
#    Y : n x m responsmatrix
#    X : n x p covariates matrix
#    R : number of bootstrap samples
#    conf : confidence level of confidence intervals
#    ests : result of GSest_multireg()


# Output:   
#   centered :  (p*m+m*m x R)-matrix of all fast/robust bootstrap recalculations
#               (recalculations are centered by original estimates)
#               (first p*m rows : regression coefficients, next m*m rows :  covariance matrix)
#               (all in vec-form, columns stacked on top of each other) 
#   
#   GSest: ((p*m + m*m) x 1) original GS estimates in vec-form
#   SE:  (p*q+q*q x 1) standard errors for GSBeta and GSSigma
#   CI.bca : (p*q+q*q x 2) 95% BCa intervals for GS-estimates
#                                                           ([lower upper])
#   CI.basic : (p*q+q*q x 2) 95% basic bootstrap intervals for GS-estimates
#                                                           ([lower upper])


 

#------------------------------------------------------------------------ 

rhobiweight <- function(x,c)
{
# Computes Tukey's biweight rho function with constant c for all values in x

hulp <- x^2/2 - x^4/(2*c^2) + x^6/(6*c^4)
rho <- hulp*(abs(x)<c) + c^2/6*(abs(x)>=c)

return(rho)
}

# --------------------------------------------------------------------

rhobiweightder1 <- function(x,c)
{
# Computes Tukey's biweight psi function met constante c voor alle waarden
# in de vector x.

hulp <- x-2*x^3/(c^2)+x^5/(c^4)
psi <- hulp*(abs(x)<c)

return(psi)
}

#-------------------------------------------------------------------------


rhobiweightder2 <- function(x,c)
{
# Computes Tukey's biweight psi function met constante c voor alle waarden
# in de vector x.

hulp <- 1-6*x^2/(c^2)+5*x^4/(c^4)
psi <- hulp*(abs(x)<c)

return(psi)
}

# --------------------------------------------------------------------

scaledpsibiweight <- function(x,c)
{
# Computes Tukey's biweight psi function with constant c for all values in x

hulp <- 1 - 2*x^2/(c^2) + x^4/(c^4)
psi <- hulp*(abs(x)<c)

return(psi)
}

#------------------------------------------------------------------------
commut <- function(p,m)
{
#computes commutation matrix
#p = no of rows
#m = no of columns (of matrix which follows Kpm)

kompm <- matrix(0,p*m,p*m)
for (k in 1:(p*m)) {
    l<-(k-1-(ceiling(k/m)-1)*m)*p+ceiling(k/m)
    kompm[k,l]<-1
}

Kpm <- kompm
return(Kpm)
}

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

#---------------------------------------------------------------------------
IRLSlocation <- function(xmat,covmat,bdp,cc)
{
xmat <- as.matrix(xmat)
n <- nrow(xmat)
p <- ncol(xmat)
neem <- sample(n,p+1)
xsub <- as.matrix(xmat[neem,])


initmu <- apply(xsub,2,mean)
initrdis <- sqrt(mahalanobis(xmat, initmu, covmat))

initobj <- mean(rhobiweight(initrdis,cc))

weights <- scaledpsibiweight(initrdis,cc)

itertest <- 0
while ((sum(weights)==0) && (itertest<500)) {
    
    neem <- sample(n,p+1)
    xsub <- as.matrix(xmat[neem,])
   
    initmu <- apply(xsub,2,mean)
    initrdis <- sqrt(mahalanobis(xmat, initmu, covmat))
    initobj <- mean(rhobiweight(initrdis,cc))
    weights <- scaledpsibiweight(initrdis,cc)
    itertest <- itertest + 1
}
if (itertest==500) stop("could not find suitable starting point for IRLS for intercept")   
   
sqrtweights <- sqrt(weights)
munieuw <- crossprod(weights, xmat) / as.vector(crossprod(sqrtweights))

rdisnieuw <-sqrt(mahalanobis(xmat,munieuw,covmat))

objnieuw <- mean(rhobiweight(rdisnieuw,cc))

iter <- 0
while (((abs(initobj/objnieuw)-1) > 10^(-15)) && (iter < 100)) {
    initobj <- objnieuw
    initmu <- munieuw
    initrdis <- rdisnieuw
    weights <- scaledpsibiweight(initrdis,cc)
    sqrtweights <- sqrt(weights)
    munieuw <- crossprod(weights, xmat) / as.vector(crossprod(sqrtweights))
    rdisnieuw <-sqrt(mahalanobis(xmat,munieuw,covmat))
    objnieuw <- mean(rhobiweight(rdisnieuw,cc))
    iter <- iter + 1
}
return(munieuw)
}

# --------------------------------------------------------------------

pvalueCI_BCA <- function(sorted, estim, infl, conf)
{
    # sorted is supposed to hold sorted and centered bootstrap estimates
    alphafunLow <- function(conf, w, a, alpha0) {
      normquan <- qnorm(1 - (1 - conf)/2)
      alphatildelow <- pnorm(w+(w-normquan)/(1-a*(w-normquan)))
      return(alphatildelow-alpha0)
    }
    Vectorize(alphafunLow)

    alphafunHigh <- function(conf, w, a, alpha0) {
      normquan <- qnorm(1 - (1 - conf)/2)
      alphatildehigh <- pnorm(w+(w+normquan)/(1-a*(w+normquan)))
      return(alphatildehigh-alpha0)
    }
    Vectorize(alphafunHigh)

    RR <- length(sorted)
    nofless <- length(sorted[sorted<=0])
    w <- qnorm(nofless/(RR+1))
    a <- 1/6 * sum(infl^3) / (sum(infl^2)^(3/2))
    alphatildelow <- alphafunLow(conf, w, a, 0)
    alphatildehigh <- alphafunHigh(conf, w, a, 0)
    indexlow <- max((RR+1)*alphatildelow,1)
    indexlow <- min(indexlow,RR)
    indexhigh <- min((RR+1)*alphatildehigh,RR)
    indexhigh <- max(indexhigh,1)
    CI <- sorted[round(indexlow)] + estim
    CI[2] <- sorted[round(indexhigh)] + estim

    # now find p-value as lowest (1-conf) for which CI includes zero:
    alpha0 <- (rank(c(0,sorted+estim))[1]-1)/(RR+1)
    if ((alpha0 == 0)||(alpha0 == RR/(RR+1))) pvalue = 0
    else {
        searchgrid = c(1e-7,1-1e-7)
        if (alphafunHigh(searchgrid[1], w, a, alpha0)<0)
            pvalue = 1-uniroot(alphafunHigh, searchgrid, w, a, alpha0)$root
        else if (alphafunLow(searchgrid[1], w, a, alpha0)>0)
            pvalue = 1-uniroot(alphafunLow, searchgrid, w, a, alpha0)$root
        else {
            pvalue = min(alpha0, 1-alpha0)*2
            warning("BCA p-value failed; simple percentile p-value is given")
        }
    }
  return(list(CI=CI, pvalue=pvalue))
}

# -----------------------------------------------------------------------------------
# -                        main function                                            -
# -----------------------------------------------------------------------------------


X <- as.matrix(X)
p<-ncol(X)
interceptdetection <- apply(X==1, 2, all)
#if (any(interceptdetection)) int=TRUE
zonderint <- (1:p)[interceptdetection==FALSE]
Xzonderint <- X[,zonderint,drop=FALSE]
X <- as.matrix(Xzonderint)

 
Y <- as.matrix(Y)
n <- nrow(Y)
m <- ncol(Y)
p <- ncol(X)

int <-nrow(ests$coefficients)> p
dimens <- p*m+m*m
dimenswith <- ifelse(int,(p+1)*m+m*m,dimens)


Im <-  diag(rep(1,m))
Ip <-  diag(rep(1,p))
 
c <- ests$c
#cc <- ests$c
b <- ests$b
bdp <-ests$method$bdp

Sigma <- ests$Sigma
if(int) {
	Betawith <- ests$coefficients
	Beta <- as.matrix(Betawith[2:(p+1),,drop=FALSE])
	Int <- t(as.matrix(Betawith[1,]))
	}
else{Beta <- ests$coefficients}	 
Sinv <- solve(Sigma)

#---------------------------------------------------------------------------
# calculate jacobian (correction matrix)                           

resmatrix <- Y-X%*%Beta

places <- t(combn(1:n,2))
ndiff <- nrow(places)
term1resvector <- as.matrix(resmatrix[places[,1],])
term2resvector <- as.matrix(resmatrix[places[,2],])
diffresmatrix <- term1resvector-term2resvector
term1X <- as.matrix(X[places[,1],])
term2X <- as.matrix(X[places[,2],])
diffX <- term1X-term2X
term1Y <- as.matrix(Y[places[,1],])
term2Y <- as.matrix(Y[places[,2],])
diffY <- term1Y-term2Y
    

divec <- sqrt(mahalanobis(diffresmatrix, rep(0,m), Sigma))
divec[divec<1e-5] <- 1e-5
udivec <- rhobiweightder1(divec,c)/divec
wdivec <- (rhobiweightder2(divec,c)*divec-rhobiweightder1(divec,c))/divec^3
zdivec <- rhobiweightder2(divec,c)
wwdivec <- rhobiweightder1(divec,c)*divec-rhobiweight(divec,c)

udipmat <- matrix(rep(udivec,p),ncol=p) * diffX
An <- crossprod(udipmat, diffX)
Vn <- crossprod(udipmat, diffY)


term1a <- matrix(0,p*p,p*m)
term1b <- matrix(0,p*m,p*m)
term2a <- matrix(0,m*m,p*m)
term2b <- matrix(0,m*m,p*m)
term2c <- matrix(0,1,p*m)
term3a <- matrix(0,p*p,m*m)
term3b <- matrix(0,p*m,m*m)
term4a <- matrix(0,m*m,m*m)
term4b <- matrix(0,1,m*m)

for (i in 1:ndiff)  {
    diffXi <- as.matrix(diffX[i,])
    diffYi <- as.matrix(diffY[i,])
    diffresi <- as.matrix(diffresmatrix[i,])
    vecdiffXidiffXi <- vecop(tcrossprod(diffXi))
    vecdiffXidiffYi <- vecop(tcrossprod(diffXi,diffYi))
    vecdiffresidiffresi <- vecop(tcrossprod(diffresi))
    wdi <- wdivec[i]
    zdi <- zdivec[i]
    udi <- udivec[i]
    tveci <- vecop(diffXi%*%t(diffresi)%*%Sinv)
    tvecSi <- vecop(Sinv%*%diffresi%*%t(diffresi)%*%Sinv)
    term1a <- term1a + wdi*tcrossprod(vecdiffXidiffXi,tveci)
    term1b <- term1b + wdi*tcrossprod(vecdiffXidiffYi,tveci)
    term2a <- term2a + wdi*tcrossprod(vecdiffresidiffresi,tveci)
    term2b <- term2b + udi*((kronecker(Im,diffresi)+kronecker(diffresi,Im))%*%(kronecker(t(diffXi),Im)%*%commut(p,m)))
    term2c <- term2c + zdi*t(tveci)
    term3a <- term3a + wdi*tcrossprod(vecdiffXidiffXi,tvecSi)
    term3b <- term3b + wdi*tcrossprod(vecdiffXidiffYi,tvecSi)
    term4a <- term4a + wdi*tcrossprod(vecdiffresidiffresi,tvecSi)
    term4b <- term4b + zdi*t(tvecSi)
}

Aninv <- solve(An)

partder1 <- (t(kronecker(Vn,Ip))%*%kronecker(t(Aninv),Aninv)%*%term1a)-(kronecker(Im,Aninv)%*%term1b)
partder2 <- -m/(b*ndiff)*(term2a+term2b)+1/(b*ndiff)*vecop(Sigma)%*%term2c
partder3 <- 1/2*(t(kronecker(Vn,Ip))%*%kronecker(t(Aninv),Aninv)%*%term3a)-1/2*(kronecker(Im,Aninv)%*%term3b)
partder4 <- -m/(2*b*ndiff)*term4a+1/(2*b*ndiff)*vecop(Sigma)%*%term4b-1/(b*ndiff)*sum(wwdivec)*diag(rep(1,m*m))

jacobian <- matrix(0,dimens,dimens)
jacobian <- cbind(rbind(partder1, partder2), rbind(partder3, partder4))


Idim <- diag(rep(1,dimens))
lincorrectmat <- solve(Idim-jacobian)

##############################################################################################

# stack Beta and Sigma into one long column vector (column by column)
vecestim <- rep(0,dimens)
vecestim[1:(p*m)] <- vecop(Beta)
vecestim[(p*m+1):dimens] <- vecop(Sigma)

# draw bootstrap samples
bootmatrix <- matrix(sample(n,R*n,replace=TRUE),ncol=R)

# now do the actual bootstrapping

bootbiasmat <- matrix(0,dimenswith,R)  
#bootbiasint <- matrix(0,m,R)  
bootsampleOK <- rep(1,R)

for (r in 1:R) {
    Yst <- as.matrix(Y[bootmatrix[,r],])
    Xst <- as.matrix(X[bootmatrix[,r],])
    term1X <- as.matrix(Xst[places[,1],])
    term2X <- as.matrix(Xst[places[,2],])
    diffXst <- term1X-term2X
    term1Y <- as.matrix(Yst[places[,1],])
    term2Y <- as.matrix(Yst[places[,2],])
    diffYst <- term1Y-term2Y
    
    diffresmatrixst <- diffYst-diffXst%*%Beta
    divecst <- sqrt(mahalanobis(diffresmatrixst, rep(0,m), Sigma))
    divecst[divecst<1e-5] <- 1e-5
    udivecst <- rhobiweightder1(divecst,c)/divecst
    wwdivecst <- rhobiweightder1(divecst,c)*divecst-rhobiweight(divecst,c)   
    
    udivecpmat <- matrix(rep(udivecst,p),ncol=p) * diffXst
    qrd <-  qr(crossprod(udivecpmat, diffXst))
    if (qrd$rank<p) {
      bootsampleOK[r] <- 0
      next
    }  
    else {
      Vnst <- solve(qrd, crossprod(udivecpmat, diffYst))
      udivecmmat <- matrix(rep(udivecst,m),ncol=m) * diffresmatrixst
      Bnst <- m * crossprod(udivecmmat, diffresmatrixst)

      vecfst <- rep(0,dimens)
      vecfst[1:(p*m)] <- vecop(Vnst) 

      Sigmast <- 1/(b*ndiff)*Bnst - 1/(b*ndiff)*sum(wwdivecst)*Sigma
      vecfst[(p*m+1):dimens] <- vecop(Sigmast)
    
      fstbias <- vecfst-vecestim
	approxest <- lincorrectmat %*% fstbias
	
if (int){	
#	approxBeta <- vecestim[1:p*m,,drop=F]+approxest[1:p*m,,drop=F]
#	approxSigma <- 
#      approxBeta <- rbind(IRLSlocation(Yst-Xst%*%reconvec(approxBeta,m),GScovariance,bdp=bdp,cc=cc),reconvec(approxBeta,m)) 
#      approxest <- c(vecop(approxBeta),approxest[-(1:p*m),,drop=F])

      approxAll <- approxest+vecestim
      approxBeta <- reconvec(approxAll[1:(p*m)],m)
      approxSigma <- reconvec(approxAll[-(1:(p*m))],m)

if (any(eigen(approxSigma,only.values = TRUE)$values<0)) approxSigma<-Sigma
	bootbiasint <-IRLSlocation(Yst-Xst%*%approxBeta,approxSigma,bdp=bdp,cc=c)- Int
	bootbiasmat[,r] <- c(bootbiasint,approxest)  
	}
else{     bootbiasmat[,r] <- approxest}  
      
    }
}

#############################################################################
if(int){vecestim=c(Int,vecestim)} 
bootindicesOK <- (1:R)[bootsampleOK==1]
ROK <- length(bootindicesOK)
nfailed <- R - ROK
if (nfailed > 0)
  warning(paste(nfailed, " out of ", R, " bootstrap samples were discarded because of too few distinct observation with positive weight"))
if (ROK>1) { 
  bootbiasmat = bootbiasmat[,bootindicesOK]
  
  #compute bootstrap estimates of variance
  GSSEs <- sqrt(apply(bootbiasmat,1,var))
  GScov <- var(t(bootbiasmat))
  # sort bootstrap recalculations for constructing intervals
  sortedGSest <- t(apply(bootbiasmat, 1, sort))

  # empirical influences for computing a in BCa intervals, based on IF(GS)
  Einf <- GSeinfs_multireg (X, Y, ests=ests)
if(int){inflE <- cbind(matrix(1,nrow=n,ncol=m),Einf$infbeta, Einf$infsigma)}
else{inflE <- cbind(Einf$infbeta, Einf$infsigma)}

  estCIbca <- matrix(0,dimenswith,2)
  estCIbasic <- matrix(0,dimenswith,2)
  pvaluebca <- rep(0,dimenswith)
  pvaluebasic <- rep(0,dimenswith)
  for (i in 1:(dimenswith)) {
    bcares <- pvalueCI_BCA(sortedGSest[i,], vecestim[i], inflE[,i], conf)
    estCIbca[i,] <- bcares$CI
    pvaluebca[i] <- bcares$pvalue
  }

  indexlow <- floor((1 - (1 - conf)/2) * ROK)
  indexhigh <- ceiling((1 - conf)/2 * ROK)
  estCIbasic[,1] <- vecestim - sortedGSest[,indexlow]
  estCIbasic[,2] <- vecestim - sortedGSest[,indexhigh]
  for (i in 1:(dimenswith)) {
    alpha.twicebeta <- (rank(c(2*vecestim[i],sortedGSest[i,]+vecestim[i]))[1]-1)/ROK
    pvaluebasic[i] <- min(alpha.twicebeta, 1-alpha.twicebeta)*2
  }
}
else
{
  warning("Too many bootstrap samples discarded; FRB is cancelled")
  bootbiasmat <- NULL
  GSSEs <- NULL
  GScov=NULL
  estCIbca <- NULL
  estCIbasic <- NULL
  pvaluebca <- NULL
  pvaluebasic <- NULL
}

return(list(centered=bootbiasmat,vecest=vecestim,SE=GSSEs,cov=GScov,
CI.bca=estCIbca,CI.basic=estCIbasic, p.bca=pvaluebca, p.basic=pvaluebasic, 
ROK=ROK))
}

##############################################################################


