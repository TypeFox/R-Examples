MMboot_twosample<-function(X, groups, R=999, ests = MMest_twosample(X,groups))
{
# robust bootstrap for two sample MM location and common shape 
# INPUT:                          
#   X : n x p data matrix
#   R : number of bootstrap samples
#   ests : result of twosampleMM
#
# OUTPUT: 
#   res$centered: (dimens x R) centered recomputations of MM and S-estimates:
#                             - first p rows: MM location center 1
#                             - second p rows: MM location center 2
#                             - next p*p rows : MM shape matrix 
#                             - next p*p rows : S covariance matrix
#                             - p rows: S location center 1
#                             - final p rows: S location center 2   
#                   (all in vec-form, columns stacked on top of each other) 
#   
#   res$MMest: (dimens x 1) original MM (and S) estimates in vec-form

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





#------------------------------------------------------------------------------
#-                          main function                                     -
#------------------------------------------------------------------------------


X<-as.matrix(X)
n<-nrow(X)
p<-ncol(X)

n1 <- sum(groups==1)
n2 <- sum(groups==2)
X1 <- X[groups==1,,drop=FALSE]   
X2 <- X[groups==2,,drop=FALSE]
dimens <- 4*p+2*p*p

mu01 <- ests$SMu1
mu02 <- ests$SMu2
Sigma0 <- ests$SSigma
MMmu1 <- ests$Mu1
MMmu2 <- ests$Mu2
MMSigma <- ests$Sigma

S0inv <- solve(Sigma0)
#auxscalesq <- det(Sigma0)^(1/p)
#auxscale <- sqrt(auxscalesq)
#MMGamma <- auxscalesq^(-1)*MMSigma


auxscale <- ests$scale
MMGamma <- ests$Gamma
MMGinv <- solve(MMGamma)

Ip<-diag(p)
c0 <- ests$c0
b <- ests$b
c1 <- ests$c1

# first fill up lower right part: the S-part 
Resi1 <- X1 - matrix(rep(mu01,n1), n1, byrow=TRUE)
divec1 <- sqrt(mahalanobis(Resi1, rep(0,p), Sigma0))
divec1[divec1<1e-5] <- 1e-5

udivec1 <- rhobiweightder1(divec1,c0)/divec1
wdivec1 <- (rhobiweightder2(divec1,c0)*divec1 - rhobiweightder1(divec1,c0))/divec1^3
zdivec1 <- rhobiweightder2(divec1,c0)
wwdivec1 <- rhobiweightder1(divec1,c0)*divec1-rhobiweight(divec1,c0)

Resi2 <- X2 - matrix(rep(mu02,n2), n2, byrow=TRUE)
divec2 <- sqrt(mahalanobis(Resi2, rep(0,p), Sigma0))
divec2[divec2<1e-5] <- 1e-5

udivec2 <- rhobiweightder1(divec2,c0)/divec2
wdivec2 <- (rhobiweightder2(divec2,c0)*divec2 - rhobiweightder1(divec2,c0))/divec2^3
zdivec2 <- rhobiweightder2(divec2,c0)
wwdivec2 <- rhobiweightder1(divec2,c0)*divec2 - rhobiweight(divec2,c0)

an1 <- sum(udivec1)
Vn1 <- crossprod(udivec1,X1)
an2 <- sum(udivec2)
Vn2 <- crossprod(udivec2,X2)

term1a <- matrix(0,1,p) 
term1b <- matrix(0,p,p)
term4a <- matrix(0,1,p*p) 
term4b <- matrix(0,p,p*p)
term6a <- matrix(0,1,p) 
term6b <- matrix(0,p,p)
term8a <- matrix(0,1,p*p) 
term8b <- matrix(0,p,p*p)
term13a <- matrix(0,p*p,p) 
term13b <- matrix(0,p*p,p) 
term13c <- matrix(0,1,p)
term14a <- matrix(0,p*p,p) 
term14b <- matrix(0,p*p,p) 
term14c <- matrix(0,1,p)
term16a <- matrix(0,p*p,p*p) 
term16b <- matrix(0,1,p*p) 
term16c <- matrix(0,p*p,p*p) 
term16d <- matrix(0,1,p*p)

for (i in 1:n1) {
    Xi <- as.matrix(X1[i,])
    resi <- as.matrix(Resi1[i,])
    vecresiresi <- vecop(resi%*%t(resi))
    wdi <- wdivec1[i]
    zdi <- zdivec1[i]
    udi <- udivec1[i]
    tveci <- crossprod(resi,S0inv)
    tvecSi <- t(vecop(S0inv%*%resi%*%t(resi)%*%S0inv))
    
    term1a <- term1a + wdi*t(resi)
    term1b <- term1b + wdi*(Xi%*%tveci)

    term4a <- term4a + wdi*tvecSi
    term4b <- term4b + wdi*(Xi%*%tvecSi)

    term13a <- term13a + wdi*(vecresiresi%*%tveci)
    term13b <- term13b + udi*(kronecker(Ip,resi)+kronecker(resi,Ip))
    term13c <- term13c + zdi*tveci
    
    term16a <- term16a + wdi*(vecresiresi%*%tvecSi)
    term16b <- term16b + zdi*tvecSi
}
for (i in 1:n2) {
    Xi <- as.matrix(X2[i,])
    resi <- as.matrix(Resi2[i,])
    vecresiresi <- vecop(resi%*%t(resi))
    wdi <- wdivec2[i]
    zdi <- zdivec2[i]
    udi <- udivec2[i]
    tveci <- crossprod(resi,S0inv)
    tvecSi <- t(vecop(S0inv%*%resi%*%t(resi)%*%S0inv))
    
    term6a <- term6a + wdi*t(resi)
    term6b <- term6b + wdi*(Xi%*%tveci)

    term8a <- term8a + wdi*tvecSi
    term8b <- term8b + wdi*(Xi%*%tvecSi)

    term14a <- term14a + wdi*(vecresiresi%*%tveci)
    term14b <- term14b + udi*(kronecker(Ip,resi)+kronecker(resi,Ip))
    term14c <- term14c + zdi*tveci
    
    term16c <- term16c + wdi*(vecresiresi%*%tvecSi)
    term16d <- term16d + zdi*tvecSi
}


partder1 <- as.matrix(t(Vn1)/an1^2)%*%term1a%*%S0inv - term1b/an1
partder2 <- matrix(0,p,p)
partder3 <- 1/2*as.matrix(t(Vn1)/an1^2)%*%term4a - 1/2*term4b/an1

partder4 <- matrix(0,p,p)
partder5 <- as.matrix(t(Vn2)/an2^2)%*%term6a%*%S0inv - term6b/an2
partder6 <- 1/2*as.matrix(t(Vn2)/an2^2)%*%term8a - 1/2*term8b/an2

partder7 <- -p/(b*n)*(term13a + term13b) + 1/(b*n)*vecop(Sigma0)%*%term13c
partder8 <- -p/(b*n)*(term14a + term14b) + 1/(b*n)*vecop(Sigma0)%*%term14c

partder9 <- -p/(2*b*n)*(term16a + term16c) + 1/(2*b*n)*vecop(Sigma0)%*%(term16b+term16d) - 1/(b*n)*(sum(wwdivec1)+sum(wwdivec2))*diag(p*p)


Part54 <- partder3
Part55 <- partder1
Part56 <- partder2

Part64 <- partder6
Part65 <- partder4
Part66 <- partder5

Part44 <- partder9
Part45 <- partder7
Part46 <- partder8



#end S-part


# now g1, g2 

Resi1 <- X1 - matrix(rep(MMmu1,n1), n1, byrow=TRUE)
divec1 <- sqrt(mahalanobis(Resi1, rep(0,p), MMGamma))
divec1[divec1<1e-5] <- 1e-5
udivec1 <- rhobiweightder1(divec1/auxscale,c1)/divec1
wdivec1 <- (rhobiweightder2(divec1/auxscale,c1)*divec1/auxscale-rhobiweightder1(divec1/auxscale,c1))/divec1^3
vdivec1 <- rhobiweightder2(divec1/auxscale,c1)

A1 <- sum(udivec1)
B1 <- crossprod(udivec1,X1)

wbig1 <- matrix(rep(sqrt(udivec1),p),ncol=p) 
wres1 <- Resi1 * wbig1  
V1 <- crossprod(wres1)


Resi2 <- X2 - matrix(rep(MMmu2,n2), n2, byrow=TRUE)
divec2 <- sqrt(mahalanobis(Resi2, rep(0,p), MMGamma))
divec2[divec2<1e-5] <- 1e-5
udivec2 <- rhobiweightder1(divec2/auxscale,c1)/divec2
wdivec2 <- (rhobiweightder2(divec2/auxscale,c1)*divec2/auxscale-rhobiweightder1(divec2/auxscale,c1))/divec2^3
vdivec2 <- rhobiweightder2(divec2/auxscale,c1)

A2 <- sum(udivec2)
B2 <- crossprod(udivec2,X2)
wbig2 <- matrix(rep(sqrt(udivec2),p),ncol=p) 
wres2 <- Resi2 * wbig2  
V2 <- crossprod(wres2)

termg1MMmua1 <- matrix(0,1,p)
termg1MMmub1 <- matrix(0,p,p)
termg1MMGammaa1 <- matrix(0,1,p*p)
termg1MMGammab1 <- matrix(0,p,p*p)
termg2MMmua1 <- matrix(0,p*p,p)
termg2MMmub1 <- matrix(0,p*p,p)
termg2MMGamma1 <- matrix(0,p*p,p*p)
termg1MMSigma0a1 <- 0
termg1MMSigma0b1 <- matrix(0,p,1)
termg2MMSigma01 <- matrix(0,p*p,1)

for (i in 1:n1) {
    Xi <- as.matrix(X1[i,])
    resi <- as.matrix(Resi1[i,])
    vecresiresi <- vecop(resi%*%t(resi))
    wdi <- wdivec1[i]
    udi <- udivec1[i]
    vdi <- vdivec1[i]
    tveci <- crossprod(resi,MMGinv)
    tvecSi <- t(vecop(MMGinv%*%resi%*%t(resi)%*%MMGinv))

    termg1MMmua1 <- termg1MMmua1+wdi*t(resi)
    termg1MMmub1 <- termg1MMmub1+wdi*(Xi%*%tveci)

    termg1MMGammaa1 <- termg1MMGammaa1+wdi*tvecSi
    termg1MMGammab1 <- termg1MMGammab1+wdi*(Xi%*%tvecSi)

    termg2MMmua1 <- termg2MMmua1+wdi*(vecresiresi%*%tveci)
    termg2MMmub1 <- termg2MMmub1+udi*(kronecker(Ip,resi)+kronecker(resi,Ip)) 

    termg2MMGamma1 <- termg2MMGamma1+wdi*(vecresiresi%*%tvecSi)
    
    termg1MMSigma0a1 <- termg1MMSigma0a1+vdi
    termg1MMSigma0b1 <- termg1MMSigma0b1+vdi*Xi

    termg2MMSigma01 <- termg2MMSigma01+vdi*vecresiresi
}    

termg1MMmua2 <- matrix(0,1,p)
termg1MMmub2 <- matrix(0,p,p)
termg1MMGammaa2 <- matrix(0,1,p*p)
termg1MMGammab2 <- matrix(0,p,p*p)
termg2MMmua2 <- matrix(0,p*p,p)
termg2MMmub2 <- matrix(0,p*p,p)
termg2MMGamma2 <- matrix(0,p*p,p*p)
termg1MMSigma0a2 <- 0
termg1MMSigma0b2 <- matrix(0,p,1)
termg2MMSigma02 <- matrix(0,p*p,1)


for (i in 1:n2) {
    Xi <- as.matrix(X2[i,])
    resi <- as.matrix(Resi2[i,])
    vecresiresi <- vecop(resi%*%t(resi))
    wdi <- wdivec2[i]
    udi <- udivec2[i]
    vdi <- vdivec2[i]
    tveci <- crossprod(resi,MMGinv)
    tvecSi <- t(vecop(MMGinv%*%resi%*%t(resi)%*%MMGinv))

    termg1MMmua2 <- termg1MMmua2+wdi*t(resi)
    termg1MMmub2 <- termg1MMmub2+wdi*(Xi%*%tveci)

    termg1MMGammaa2 <- termg1MMGammaa2+wdi*tvecSi
    termg1MMGammab2 <- termg1MMGammab2+wdi*(Xi%*%tvecSi)

    termg2MMmua2 <- termg2MMmua2+wdi*(vecresiresi%*%tveci)
    termg2MMmub2 <- termg2MMmub2+udi*(kronecker(Ip,resi)+kronecker(resi,Ip)) 

    termg2MMGamma2 <- termg2MMGamma2+wdi*(vecresiresi%*%tvecSi)
    
    termg1MMSigma0a2 <- termg1MMSigma0a2+vdi
    termg1MMSigma0b2 <- termg1MMSigma0b2+vdi*Xi

    termg2MMSigma02 <- termg2MMSigma02+vdi*vecresiresi
}    

Part11 <- as.matrix(t(B1)/A1^2)%*%termg1MMmua1%*%MMGinv-termg1MMmub1/A1
Part22 <- as.matrix(t(B2)/A2^2)%*%termg1MMmua2%*%MMGinv-termg1MMmub2/A2
Part13 <- 1/2*as.matrix(t(B1)/A1^2)%*%termg1MMGammaa1-1/2*termg1MMGammab1/A1
Part23 <- 1/2*as.matrix(t(B2)/A2^2)%*%termg1MMGammaa2-1/2*termg1MMGammab2/A2
Part14 <- -1/2/p/auxscale*(as.matrix(t(B1)/A1^2)%*%termg1MMSigma0a1-termg1MMSigma0b1/A1)%*%t(vecop(t(S0inv)))
Part24 <- -1/2/p/auxscale*(as.matrix(t(B2)/A2^2)%*%termg1MMSigma0a2-termg1MMSigma0b2/A2)%*%t(vecop(t(S0inv)))
Part31 <- -det((V1+V2))^(-1/p)*(diag(p*p)-1/p*vecop(V1+V2)%*%t(vecop(t(solve(V1+V2)))))%*%(termg2MMmua1+termg2MMmub1)
Part32 <- -det((V1+V2))^(-1/p)*(diag(p*p)-1/p*vecop(V1+V2)%*%t(vecop(t(solve(V1+V2)))))%*%(termg2MMmua2+termg2MMmub2)
Part33 <- -1/2*det(V1+V2)^(-1/p)*(diag(p*p)-1/p*vecop(V1+V2)%*%t(vecop(t(solve(V1+V2)))))%*%(termg2MMGamma1+termg2MMGamma2)
Part34 <- -1/2/p/auxscale*det(V1+V2)^(-1/p)*(diag(p*p)-1/p*vecop(V1+V2)%*%t(vecop(t(solve(V1+V2)))))%*%(termg2MMSigma01+termg2MMSigma02)%*%t(vecop(t(S0inv)))

Part12 <- matrix(0,p,p)
Part15 <- matrix(0,p,p)
Part16 <- matrix(0,p,p)
Part21 <- matrix(0,p,p)
Part25 <- matrix(0,p,p)
Part26 <- matrix(0,p,p)
Part35 <- matrix(0,p*p,p)
Part36 <- matrix(0,p*p,p)
Part41 <- matrix(0,p*p,p)
Part42 <- matrix(0,p*p,p)
Part43 <- matrix(0,p*p,p*p)
Part51 <- matrix(0,p,p)
Part52 <- matrix(0,p,p)
Part53 <- matrix(0,p,p*p)
Part61 <- matrix(0,p,p)
Part62 <- matrix(0,p,p)
Part63 <- matrix(0,p,p*p)

col1 <- rbind(Part11, Part21, Part31, Part41,Part51,Part61)
col2 <- rbind(Part12, Part22, Part32, Part42,Part52,Part62)
col3 <- rbind(Part13, Part23, Part33, Part43,Part53,Part63)
col4 <- rbind(Part14, Part24, Part34, Part44,Part54,Part64)
col5 <- rbind(Part15, Part25, Part35, Part45,Part55,Part65)
col6 <- rbind(Part16, Part26, Part36, Part46,Part56,Part66)

jacobian <- cbind(col1, col2, col3, col4,col5,col6)


Idim <- diag(dimens)
lincorrectmat <- solve(Idim-jacobian)


# now do the actual bootstrapping

# put all parameters into one large vector:
vecestim <- rep(0,dimens)
vecestim[1:p] <- vecop(MMmu1)
vecestim[(p+1):(2*p)] <- vecop(MMmu2)
vecestim[(2*p+1):(2*p+p*p)] <- vecop(MMGamma)
vecestim[(2*p+p*p+1):(2*p+2*p*p)] <- vecop(Sigma0)
vecestim[(2*p+2*p*p+1):(3*p+2*p*p)] <- vecop(mu01)
vecestim[(3*p+2*p*p+1):(4*p+2*p*p)]<- vecop(mu02)

# to draw bootstrap samples
#set.seed(2)
bootmatrix1 <- matrix(sample(n1,R*n1,replace=TRUE),ncol=R)
bootmatrix2 <- matrix(sample(n2,R*n2,replace=TRUE),ncol=R)



bootbiasmat <- matrix(0,dimens,R)  
bootbiaszc <- matrix(0,dimens,R) 

for (r in 1:R) {
    X1st <- X1[bootmatrix1[,r],,drop=FALSE]
    Resi1st <- X1st - matrix(rep(MMmu1,n1), n1, byrow=TRUE)
    divec1st <- sqrt(mahalanobis(Resi1st, rep(0,p), MMSigma))
    divec1st[divec1st<1e-5] <- 1e-5
    udivec1st <- rhobiweightder1(divec1st,c1)/divec1st
    
    restildematrix1st <- X1st-matrix(rep(mu01,n1), n1, byrow=TRUE)
    ditildevec1st <- sqrt(mahalanobis(restildematrix1st,rep(0,p),Sigma0))
    ditildevec1st[ditildevec1st<1e-5] <- 1e-5
    uditildevec1st <- rhobiweightder1(ditildevec1st,c0)/ditildevec1st
    wwditildevec1st=rhobiweightder1(ditildevec1st,c0)*ditildevec1st-rhobiweight(ditildevec1st,c0)   
    
    B1st=(1/sum(udivec1st))*crossprod(udivec1st,X1st)

    sqrtwu1 <- sqrt(udivec1st)
    wbig1 <- matrix(rep(sqrtwu1,p),ncol=p) 
    wres1 <- Resi1st * wbig1  
    G1st <- crossprod(wres1)
    
    sqrtwu1 <- sqrt(uditildevec1st)
    wbig1 <- matrix(rep(sqrtwu1,p),ncol=p) 
    wres1 <- restildematrix1st * wbig1  
    V01stterm1 <- 1/(b*n)*p*crossprod(wres1)
    V01stterm2 <- 1/(b*n)*sum(wwditildevec1st)*Sigma0
    V01st <- V01stterm1-V01stterm2
 
    B01st <- (1/sum(uditildevec1st))*crossprod(uditildevec1st,X1st)

    X2st <- X2[bootmatrix2[,r],,drop=FALSE]
    Resi2st <- X2st - matrix(rep(MMmu2,n2), n2, byrow=TRUE)
    divec2st <- sqrt(mahalanobis(Resi2st, rep(0,p), MMSigma))
    divec2st[divec2st<1e-5] <- 1e-5
    udivec2st <- rhobiweightder1(divec2st,c1)/divec2st
    
    restildematrix2st <- X2st-matrix(rep(mu02,n2), n2, byrow=TRUE)
    ditildevec2st <- sqrt(mahalanobis(restildematrix2st,rep(0,p),Sigma0))
    ditildevec2st[ditildevec2st<1e-5] <- 1e-5
    uditildevec2st <- rhobiweightder1(ditildevec2st,c0)/ditildevec2st
    wwditildevec2st <- rhobiweightder1(ditildevec2st,c0)*ditildevec2st-rhobiweight(ditildevec2st,c0)   
    
    B2st<-(1/sum(udivec2st))*crossprod(udivec2st,X2st)
    sqrtwu2 <- sqrt(udivec2st)
    wbig2 <- matrix(rep(sqrtwu2,p),ncol=p) 
    wres2 <- Resi2st * wbig2  
    G2st <- crossprod(wres2)
    
    sqrtwu2 <- sqrt(uditildevec2st)
    wbig2 <- matrix(rep(sqrtwu2,p),ncol=p) 
    wres2 <- restildematrix2st * wbig2  
    V02stterm1 <- 1/(b*n)*p*crossprod(wres2)
    V02stterm2 <- 1/(b*n)*sum(wwditildevec2st)*Sigma0
    V02st <- V02stterm1-V02stterm2
    
    B02st <- (1/sum(uditildevec2st))*crossprod(uditildevec2st,X2st)
    
    vecfst <- rep(0,dimens)
    vecfst[1:p] <- vecop(B1st)
    vecfst[(p+1):(2*p)] <- vecop(B2st)
    Gammast <- det((G1st+G2st))^(-1/p)*(G1st+G2st)
    vecfst[(2*p+1):(2*p+p*p)] <- vecop(Gammast)
    
    Sigmast <- V01st+V02st
    vecfst[(2*p+p*p+1):(2*p+2*p*p)] <- vecop(Sigmast)
    vecfst[(2*p+2*p*p+1):(3*p+2*p*p)] <- vecop(B01st)
    vecfst[(3*p+2*p*p+1):(4*p+2*p*p)] <- vecop(B02st)
    
    fstbias <- vecfst - vecestim
    
    bootbiaszc[,r] <- fstbias
    
    # apply linear correction:
    bootbiasmat[,r] <- lincorrectmat %*% fstbias 
}    
    
return(list(centered=bootbiasmat,MMest=vecestim))
}

