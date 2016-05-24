Sboot_multireg<- function(X, Y, R=999, conf=0.95, ests=Sest_multireg(X, Y))
{
# robust bootstrap for multivariate S-regression
# INPUT:
#   Y : n x q response matrix
#   X : n x p covariates matrix (ones(n,1) for location/shape estimation)
#   R : number of bootstrap samples
#   conf : confidence level for bootstrap intervals
#   ests : result of Sest_multireg()
#
# OUTPUT: 
#   res$centered: ((p*q + q*q) x R) centered recomputations S-estimates:
#                             - first p*q rows: regression coeffients      
#                             - next q*q rows : covariance matrix
#                   (all in vec-form, columns stacked on top of each other) 
#   
#   res$vecest: ((p*q + q*q) x 1) original S estimates in vec-form
#   res$SE: ((p*q + q*q) x 1) bootstrap standard errors for elements in res$Sest
#   res$CI.bca: ((p*q + q*q) x 2) BCa confidence limits for elements in res$Sest
#   res$CI.basic: ((p*q + q*q) x 2) basic bootstrap confidence limits for elements in res$Sest
# --------------------------------------------------------------------

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
# Computes Tukey's biweight psi function with constant c for all values in x

hulp <- x - 2*x^3/(c^2) + x^5/(c^4)
rho <- hulp*(abs(x)<c)

return(rho)
}

# --------------------------------------------------------------------

rhobiweightder2 <- function(x,c)
{
# Computes derivative of Tukey's biweight psi function with constant c for all values in x

hulp <- 1 - 6*x^2/(c^2) + 5*x^4/(c^4)
rho <- hulp*(abs(x)<c)

return(rho)
}

# --------------------------------------------------------------------

commut <- function(p,m) {

# computes commutation matrix
# p = no. of rows
# m = no. of columns (of matrix which follows the commut matrix)

kompm <- matrix(0,p*m,p*m)
for (k in 1:(p*m)) {
    l <- (k - 1 - (ceiling(k/m)-1) * m) * p + ceiling(k/m)
    kompm[k,l] <- 1
}

return(kompm)

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
# --------------------------------------------------------------------
# -                         main function                            -
# --------------------------------------------------------------------

Y <- as.matrix(Y)
X <- as.matrix(X)
n <- nrow(X)
if (nrow(ests$coefficients)>ncol(X)) X <- cbind(rep(1,n),X)
p <- ncol(X)
q <- ncol(Y)
dimens <- p*q + q*q

Iq <- diag(rep(1,q))
Ip <- diag(rep(1,p))

c <- ests$c
b <- ests$b
Sigma0 <- ests$Sigma
Beta0 <- ests$coefficients

Sinv <- solve(Sigma0)

###############################################################
###                 calculate jacobian                      ###
###############################################################

#           p*q          q*q        
#       -------------------------
#       |           |            |
#  p*q  | g1_Beta   | g1_Sigma   |
#       |           |            |
#       -------------------------
#       |           |            |
#  q*q  | g2_Beta   | g2_Sigma   |
#       |           |            |  
#       -------------------------

resmatrix <- Y - X %*% Beta0
divec <- sqrt(mahalanobis(resmatrix, rep(0,q), Sigma0))
divec[divec < 1e-5] <- 1e-5
udivec <- rhobiweightder1(divec,c)/divec
wdivec <- (rhobiweightder2(divec,c)*divec - rhobiweightder1(divec,c))/divec^3   
zdivec <- rhobiweightder2(divec,c)
wwdivec <- rhobiweightder1(divec,c)*divec - rhobiweight(divec,c)

#A <- t(X) %*% (matrix(rep(udivec,p),ncol=p) * X)
#B <- t(X) %*% (matrix(rep(udivec,q),ncol=q) * Y)
uX <- matrix(rep(udivec,p),ncol=p) * X
A <- crossprod(uX, X)
B <- crossprod(uX, Y)

term1a <- matrix(0, p*p,p*q)
term1b <- matrix(0,p*q,p*q)
term2a <- matrix(0,q*q,p*q)
term2b <- matrix(0,q*q,p*q)
term2c <- matrix(0,1,p*q)
term3a <- matrix(0,p*p,q*q)
term3b <- matrix(0,p*q,q*q)
term4a <- matrix(0,q*q,q*q)
term4b <- matrix(0,1,q*q)

for (i in 1:n) {
    Xi <- as.matrix(X[i,])
    Yi <- as.matrix(Y[i,])
    resi <- as.matrix(resmatrix[i,])
#    vecXiXi <- vecop(Xi %*% t(Xi))
#    vecXiYi <- vecop(Xi %*% t(Yi))
    vecXiXi <- vecop(tcrossprod(Xi))
    vecXiYi <- vecop(tcrossprod(Xi,Yi))
    vecresiresi <- vecop(tcrossprod(resi))
    wdi <- wdivec[i]
    zdi <- zdivec[i]
    udi <- udivec[i]
    veci <- vecop(Xi %*% t(resi) %*% Sinv)
    vecSi <- vecop(Sinv %*% resi %*% t(resi) %*% Sinv)

    term1a <- term1a + wdi * tcrossprod(vecXiXi, veci)
    term1b <- term1b + wdi * tcrossprod(vecXiYi, veci)

    term2a <- term2a + wdi * tcrossprod(vecresiresi, veci)
    term2b <- term2b + udi * ((kronecker(Iq,resi) + kronecker(resi,Iq)) %*% (kronecker(t(Xi),Iq) %*% commut(p,q))) 
    term2c <- term2c + zdi * t(veci)
    
    term3a <- term3a + wdi * tcrossprod(vecXiXi, vecSi)
    term3b <- term3b + wdi * tcrossprod(vecXiYi, vecSi)

    term4a <- term4a + wdi * tcrossprod(vecresiresi, vecSi)
    term4b <- term4b + zdi * t(vecSi)
}

Ainv <- solve(A)

partder1 <- (t(kronecker(B,Ip)) %*% kronecker(t(Ainv),Ainv) %*% term1a) - (kronecker(Iq,Ainv) %*% term1b)   
partder2 <- -q/(b*n) * (term2a + term2b) + 1/(b*n) * vecop(Sigma0) %*% term2c
partder3 <- 1/2 * (t(kronecker(B,Ip)) %*% kronecker(t(Ainv),Ainv) %*% term3a) - 1/2 * (kronecker(Iq,Ainv) %*% term3b)  
partder4 <- -q/(2*b*n) * term4a + 1/(2*b*n) * vecop(Sigma0) %*% term4b - 1/(b*n) * sum(wwdivec) * diag(rep(1,q*q))

jacobian <- cbind(rbind(partder1, partder2), rbind(partder3, partder4))

Idim <- diag(rep(1,dimens))
lincorrectmat <- solve(Idim-jacobian)

######################################################################

# put all estimates (coefs and covariances) in one column 
vecestim <- rep(0,dimens)
vecestim[1:(p*q)] <- vecop(Beta0)
vecestim[(p*q+1):dimens] <- vecop(Sigma0)

# to draw bootstrap samples 
#set.seed(2)
bootmatrix <- matrix(sample(n,R*n,replace=TRUE),ncol=R)

bootbiasmat <- matrix(0,dimens,R)  
bootsampleOK <- rep(1,R)

for (r in 1:R) {
    bootind <- bootmatrix[,r]
    Yst <- Y[bootind,]
    Xst <- X[bootind,]
    resmatrixst <- resmatrix[bootind,]
    udivecst <- udivec[bootind]
    wwdivecst <- wwdivec[bootind]
    
    uXst <- matrix(rep(udivecst,p),ncol=p) * Xst
    qrd <-  qr(crossprod(uXst, Xst))
    if (qrd$rank<p) {
      bootsampleOK[r] <- 0
      next
    }  
    else {
      Bst <- solve(qrd, crossprod(uXst, Yst))
      uresst <- matrix(rep(udivecst,q),ncol=q) * resmatrixst
      Vst_term1 <- 1/(b*n) * q * crossprod(uresst, resmatrixst)
      Vst_term2 <- 1/(b*n)* sum(wwdivecst) * Sigma0
      Vst <- Vst_term1 - Vst_term2
    
      # list uncorrected bootstrap recomputations
      vecfst <- rep(0,dimens)
      vecfst[1:(p*q)] <- vecop(Bst)
      vecfst[(p*q+1):dimens] <- vecop(Vst)
    
      # compute centered, corrected fast bootstrap estimates
      fstbias <- vecfst - vecestim
      bootbiasmat[,r] <- lincorrectmat %*% fstbias  
    }
}

#############################################################################

bootindicesOK <- (1:R)[bootsampleOK==1]
ROK <- length(bootindicesOK)
nfailed <- R - ROK
if (nfailed > 0)
  warning(paste(nfailed, " out of ", R, " bootstrap samples were discarded because of too few distinct observation with positive weight"))
if (ROK>1) { 
  bootbiasmat = bootbiasmat[,bootindicesOK]
  
  # compute bootstrap estimates of standard error
  SSEs <- sqrt(apply(bootbiasmat, 1, var))
  Scov <- var(t(bootbiasmat))
  # sort bootstrap recalculations for constructing intervals
  sortedSest <- t(apply(bootbiasmat, 1, sort))
  
  # empirical inlfuences for computing a in BCa intervals, based on IF(S)
  Einf <- Seinfs_multireg(X, Y, ests=ests)
  inflE <- cbind(Einf$BetaS, Einf$covS)
  
  estCIbca <- matrix(0,dimens,2)
  estCIbasic <- matrix(0,dimens,2)
  pvaluebca <- rep(0,dimens)
  pvaluebasic <- rep(0,dimens)
  for (i in 1:(dimens)) {
    bcares <- pvalueCI_BCA(sortedSest[i,], vecestim[i], inflE[,i], conf)
    estCIbca[i,] <- bcares$CI
    pvaluebca[i] <- bcares$pvalue
  }
  
  indexlow <- floor((1 - (1 - conf)/2) * ROK)
  indexhigh <- ceiling((1 - conf)/2 * ROK)
  estCIbasic[,1] <- vecestim - sortedSest[,indexlow]
  estCIbasic[,2] <- vecestim - sortedSest[,indexhigh]
  for (i in 1:(dimens)) {
    alpha.twicebeta <- (rank(c(2*vecestim[i],sortedSest[i,]+vecestim[i]))[1]-1)/ROK
    pvaluebasic[i] <- min(alpha.twicebeta, 1-alpha.twicebeta)*2
  }
}
else
{
  warning("Too many bootstrap samples discarded; FRB is cancelled")
  bootbiasmat <- NULL
  SSEs <- NULL
  Scov=NULL
  estCIbca <- NULL
  estCIbasic <- NULL
  pvaluebca <- NULL
  pvaluebasic <- NULL
  }
#############################################################################

return(list(centered=bootbiasmat, vecest=vecestim, SE=SSEs, cov=Scov, 
CI.bca=estCIbca, p.bca=pvaluebca, p.basic=pvaluebasic, CI.basic=estCIbasic, 
ROK=ROK))

}

