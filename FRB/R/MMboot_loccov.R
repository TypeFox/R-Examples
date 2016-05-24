MMboot_loccov <- function(Y, R=999, ests = MMest_loccov(Y))
{
# robust bootstrap for multivariate MM location/shape estimation 
# INPUT:
#   Y : n x q data matrix
#   R : number of bootstrap samples
#   ests : result of multiMM_location
#
# OUTPUT: 
#   res$centered: (2*dimens x R) centered recomputations of MM and S-estimates:
#                             - first q rows: MM location 
#                             - next q*q rows : MM shape matrix 
#                             - next q*q rows : S covariance matrix
#                             - final q rows: S location 
#                   (all in vec-form, columns stacked on top of each other) 
#   
#   res$MMest: (2*dimens x 1) original MM (and S) estimates in vec-form

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
# -                         main function                            -
# --------------------------------------------------------------------

Y <- as.matrix(Y)
n <- nrow(Y)
q <- ncol(Y)
dimens <- q + q*q

Iq <- diag(rep(1,q))

c0 <- ests$c0
b <- ests$b
c1 <- ests$c1

MMMu <- ests$Mu
MMSigma <- ests$Sigma
Sigma0 <- ests$SSigma
Mu0 <- ests$SMu

MMSinv <- solve(MMSigma)
S0inv <- solve(Sigma0)
auxscalesq <- det(Sigma0)^(1/q)
auxscale <- sqrt(auxscalesq)
MMGamma <- auxscalesq^(-1)*MMSigma
MMGinv <- solve(MMGamma)

##########################################################################
###                        calculate jacobian                          ###
##########################################################################

#             q          q*q          q*q            q  
#       -----------------------------------------------------
#       |           |            |            |             |
#    q  | g1_MMBeta | g1_MMGamma | g1_Sigma_0 |      0      |
#       |           |            |            |             |
#       -----------------------------------------------------
#       |           |            |            |             |
#  q*q  | g2_MMBeta | g2_MMGamma | g2_Sigma_0 |      0      |
#       |           |            |            |             |
#       -----------------------------------------------------
#       |           |            |            |             |
#  q*q  |     0     |      0     | g3_Sigma_0 |  g3_Beta_0  |
#       |           |            |            |             |
#       -----------------------------------------------------
#       |           |            |            |             |
#    q  |     0     |      0     | g4_Sigma_0 |  g4_Beta_0  |
#       |           |            |            |             |
#       -----------------------------------------------------
#             1            2            3            4 

# first fill up lower right part: the S-part : g3, g4

restildematrix <- Y - matrix(rep(Mu0,n), n, byrow=TRUE)
ditildevec <- sqrt(mahalanobis(restildematrix, rep(0,q), Sigma0))
ditildevec[ditildevec < 1e-5] <- 1e-5
uditildevec <- rhobiweightder1(ditildevec,c0)/ditildevec
wditildevec <- (rhobiweightder2(ditildevec,c0)*ditildevec - rhobiweightder1(ditildevec,c0))/ditildevec^3   
zditildevec <- rhobiweightder2(ditildevec,c0)
wwditildevec <- rhobiweightder1(ditildevec,c0)*ditildevec - rhobiweight(ditildevec,c0)

atilde <- sum(uditildevec)
btilde <- crossprod(uditildevec, Y)

term1a <- matrix(0,1,q)
term1b <- matrix(0,q,q)
term2a <- matrix(0,q*q,q)
term2b <- matrix(0,q*q,q)
term2c <- matrix(0,1,q)
term3a <- matrix(0,1,q*q)
term3b <- matrix(0,q,q*q)
term4a <- matrix(0,q*q,q*q)
term4b <- matrix(0,1,q*q)

for (i in 1:n) {
    Yi <- as.matrix(Y[i,])
    resi <- as.matrix(restildematrix[i,])
    vecresiresi <- vecop(resi %*% t(resi))
    wdi <- wditildevec[i]
    zdi <- zditildevec[i]
    udi <- uditildevec[i]
    tveci <- crossprod(resi, S0inv)
    tvecSi <- t(vecop(S0inv %*% resi %*% t(resi) %*% S0inv))

    term1a <- term1a + wdi * t(resi)
    term1b <- term1b + wdi * (Yi %*% tveci)

    term2a <- term2a + wdi * (vecresiresi %*% tveci)
    term2b <- term2b + udi * (kronecker(Iq,resi) + kronecker(resi,Iq)) 
    term2c <- term2c + zdi * tveci
    
    term3a <- term3a + wdi * tvecSi
    term3b <- term3b + wdi * (Yi %*% tvecSi)

    term4a <- term4a + wdi * (vecresiresi %*% tvecSi)
    term4b <- term4b + zdi * tvecSi
}

partder1 <- as.matrix(t(btilde) / atilde^2)  %*% term1a %*% S0inv - term1b / atilde   
partder2 <- -q/(b*n) * (term2a + term2b) + 1/(b*n) * vecop(Sigma0) %*% term2c
partder3 <- 1/2 * as.matrix(t(btilde) / atilde^2) %*% term3a - 1/2 * term3b / atilde 
partder4 <- -q/(2*b*n) * term4a + 1/(2*b*n) * vecop(Sigma0) %*% term4b - 1/(b*n) * sum(wwditildevec) * diag(rep(1,q*q))

Part33 <- partder4
Part34 <- partder2
Part43 <- partder3
Part44 <- partder1

# end S-part

# now g1, g2 

resmatrix <- Y - matrix(rep(MMMu,n), n, byrow=TRUE)
divec <- sqrt(mahalanobis(resmatrix, rep(0,q), MMGamma))
divec[divec < 1e-5] <- 1e-5
udivec <- rhobiweightder1(divec/auxscale,c1)/divec
wdivec <- (rhobiweightder2(divec/auxscale,c1)*divec/auxscale - rhobiweightder1(divec/auxscale,c1))/divec^3 
vdivec <- rhobiweightder2(divec/auxscale,c1)

a <- sum(udivec)
B <- crossprod(udivec, Y)
wbig <- matrix(rep(sqrt(udivec),q),ncol=q) 
wres <- resmatrix * wbig  
V <- crossprod(wres)

termg1MMBetaa <- matrix(0,1,q);
termg1MMBetab <- matrix(0,q,q);
termg1MMGammaa <- matrix(0,1,q*q);
termg1MMGammab <- matrix(0,q,q*q);
termg2MMBetaa <- matrix(0,q*q,q);
termg2MMBetab <- matrix(0,q*q,q);
termg2MMGamma <- matrix(0,q*q,q*q);
termg1MMSigma0a <- 0;
termg1MMSigma0b <- matrix(0,q,1);
termg2MMSigma0 <- matrix(0,q*q,1);

for (i in 1:n) {
    Yi <- as.matrix(Y[i,])
    resi <- as.matrix(resmatrix[i,])
    vecresiresi <- vecop(resi %*% t(resi))
    wdi <- wdivec[i]
    udi <- udivec[i]
    vdi <- vdivec[i]
    tveci <- crossprod(resi, MMGinv)
    tvecSi <- t(vecop(MMGinv %*% resi %*% t(resi) %*% MMGinv))

    termg1MMBetaa <- termg1MMBetaa + wdi * t(resi)
    termg1MMBetab <- termg1MMBetab + wdi * (Yi %*% tveci)

    termg1MMGammaa <- termg1MMGammaa + wdi * tvecSi
    termg1MMGammab <- termg1MMGammab + wdi * (Yi %*% tvecSi)

    termg2MMBetaa <- termg2MMBetaa + wdi * (vecresiresi %*% tveci)
    termg2MMBetab <- termg2MMBetab + udi * (kronecker(Iq,resi) + kronecker(resi,Iq)) 

    termg2MMGamma <- termg2MMGamma + wdi * (vecresiresi %*% tvecSi)
    
    termg1MMSigma0a <- termg1MMSigma0a + vdi
    termg1MMSigma0b <- termg1MMSigma0b + vdi * Yi

    termg2MMSigma0 <- termg2MMSigma0 + vdi * vecresiresi
}

Part11 <- as.matrix(t(B) / a^2) %*% termg1MMBetaa %*% MMGinv - termg1MMBetab / a
Part12 <- 1/2 * as.matrix(t(B) / a^2) %*% termg1MMGammaa - 1/2 * termg1MMGammab / a
Part21 <- -det(V)^(-1/q) * (diag(rep(1,q*q)) - 1/q * vecop(V) %*% t(vecop(t(solve(V))))) %*% (termg2MMBetaa + termg2MMBetab)
Part22 <- -1/2 * det(V)^(-1/q) * (diag(rep(1,q*q)) - 1/q * vecop(V) %*% t(vecop(t(solve(V))))) %*% termg2MMGamma
Part13 <- -1/2/q/auxscale * (as.matrix(t(B) / a^2) %*% termg1MMSigma0a - termg1MMSigma0b / a) %*% t(vecop(t(S0inv)))   
Part23 <- -1/2/q/auxscale * det(V)^(-1/q) * (diag(rep(1,q*q)) - 1/q * vecop(V) %*% t(vecop(t(solve(V))))) %*% termg2MMSigma0 %*% t(vecop(t(S0inv)))

Part14 <- matrix(0,q,q)
Part24 <- matrix(0,q*q,q)
Part31 <- matrix(0,q*q,q)
Part32 <- matrix(0,q*q,q*q)
Part41 <- matrix(0,q,q)
Part42 <- matrix(0,q,q*q)

col1 <- rbind(Part11, Part21, Part31, Part41)
col2 <- rbind(Part12, Part22, Part32, Part42)
col3 <- rbind(Part13, Part23, Part33, Part43)
col4 <- rbind(Part14, Part24, Part34, Part44)
jacobian <- cbind(col1, col2, col3, col4)

Idim <- diag(rep(1,dimens*2))
lincorrectmat <- solve(Idim-jacobian)

######################################################################

# put all estimates (coefs and covariances) in one column 
vecestim <- rep(0,dimens*2)
vecestim[1:q] <- vecop(MMMu)
vecestim[(q+1):dimens] <- vecop(MMGamma)
vecestim[(dimens+1):(dimens+(q*q))] <- vecop(Sigma0)
vecestim[(dimens+q*q+1):(dimens*2)] <- vecop(Mu0)

# to draw bootstrap samples
#set.seed(2)
bootmatrix <- matrix(sample(n,R*n,replace=TRUE),ncol=R)

bootbiasmat <- matrix(0,dimens*2,R)  

for (r in 1:R) {
    Yst <- Y[bootmatrix[,r],]
    resmatrixst <- Yst - matrix(rep(MMMu,n), n, byrow=TRUE)
    divecst <- sqrt(mahalanobis(resmatrixst, rep(0,q), MMSigma))
    divecst[divecst<1e-5] <- 1e-5
    udivecst <- rhobiweightder1(divecst,c1)/divecst
    
    restildematrixst <- Yst - matrix(rep(Mu0,n), n, byrow=TRUE) 
    ditildevecst<- sqrt(mahalanobis(restildematrixst, rep(0,q), Sigma0))
    ditildevecst[ditildevecst<1e-5] <- 1e-5
    uditildevecst <- rhobiweightder1(ditildevecst,c0)/ditildevecst
    wwditildevecst <- rhobiweightder1(ditildevecst,c0)*ditildevecst - rhobiweight(ditildevecst,c0)     
    
    sqrtwu <- sqrt(udivecst)
    Bst <- crossprod(udivecst, Yst) / as.vector(crossprod(sqrtwu))
    wbig <- matrix(rep(sqrtwu,q),ncol=q) 
    wres <- resmatrixst * wbig  
    Gst <- crossprod(wres)
    Gst <- (determinant(Gst,logarithm=FALSE)$modulus)^(-1/q) * Gst  
    Vst <- auxscalesq * Gst
    
    sqrtwu <- sqrt(uditildevecst)
    B0st <- crossprod(uditildevecst, Yst) / as.vector(crossprod(sqrtwu))
    wbig <- matrix(rep(sqrtwu,q),ncol=q) 
    wres <- restildematrixst * wbig  
    V0st_term1 <- 1/(b*n)*q * crossprod(wres)
    V0st_term2 <- 1/(b*n)* sum(wwditildevecst) * Sigma0
    V0st <- V0st_term1 - V0st_term2
    
    # list uncorrected bootstrap recomputations
    vecfst <- rep(0,dimens*2)
    vecfst[1:q] <- vecop(Bst)
    vecfst[(q+1):dimens] <- vecop(Gst)
    vecfst[(dimens+1):(dimens+q*q)] <- vecop(V0st)
    vecfst[(dimens+q*q+1):(dimens*2)] <- vecop(B0st)
    
    # compute centered, corrected fast bootstrap estimates
    fstbias <- vecfst - vecestim
    bootbiasmat[,r] <- lincorrectmat %*% fstbias  
 
}

##########################################################################

return(list(centered=bootbiasmat, MMest=vecestim))

}
