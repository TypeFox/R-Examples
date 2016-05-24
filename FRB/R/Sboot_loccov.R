Sboot_loccov <- function(Y, R=999, ests=Sest_loccov(Y))
{
# robust bootstrap for multivariate S location covariance estimation 
# INPUT:
#   Y : n x q data matrix
#   R : number of bootstrap samples
#   ests : result of Sest_loccov
#
# OUTPUT: 
#   res$centered: ((q + q*q) x R) centered recomputations of S-estimates:
#                             - first q rows: S location 
#                             - next q*q rows : S covariance matrix 
#                   (in vec-form, columns stacked on top of each other) 
#   
#   res$Sest: ((q + q*q) x 1) original S estimates in vec-form

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

c0 <- ests$c
b <- ests$b

Sigma0 <- ests$Sigma
Mu0 <- ests$Mu
S0inv <- solve(Sigma0)

##########################################################################
###                        calculate jacobian                          ###
##########################################################################

#            q          q*q        
#       -------------------------
#       |           |            |
#    q  | g1_Mu     | g1_Sigma   |
#       |           |            |
#       -------------------------
#       |           |            |
#  q*q  | g2_Mu     | g2_Sigma   |
#       |           |            |  
#       -------------------------

resmatrix <- Y - matrix(rep(Mu0,n), n, byrow=TRUE)
divec <- sqrt(mahalanobis(resmatrix, rep(0,q), Sigma0))
divec[divec < 1e-5] <- 1e-5
udivec <- rhobiweightder1(divec,c0)/divec
wdivec <- (rhobiweightder2(divec,c0)*divec - rhobiweightder1(divec,c0))/divec^3   
zdivec <- rhobiweightder2(divec,c0)
wwdivec <- rhobiweightder1(divec,c0)*divec - rhobiweight(divec,c0)

a <- sum(udivec)
B <- crossprod(udivec, Y)

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
    resi <- as.matrix(resmatrix[i,])
    vecresiresi <- vecop(resi %*% t(resi))
    wdi <- wdivec[i]
    zdi <- zdivec[i]
    udi <- udivec[i]
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

partder1 <- as.matrix(t(B) / a^2)  %*% term1a %*% S0inv - term1b / a   
partder2 <- -q/(b*n) * (term2a + term2b) + 1/(b*n) * vecop(Sigma0) %*% term2c
partder3 <- 1/2 * as.matrix(t(B) / a^2) %*% term3a - 1/2 * term3b / a 
partder4 <- -q/(2*b*n) * term4a + 1/(2*b*n) * vecop(Sigma0) %*% term4b - 1/(b*n) * sum(wwdivec) * diag(rep(1,q*q))

jacobian <- cbind(rbind(partder1, partder2), rbind(partder3, partder4))

Idim <- diag(rep(1,dimens))
lincorrectmat <- solve(Idim-jacobian)

######################################################################

# put all estimates (coefs and covariances) in one column 
vecestim <- rep(0,dimens)
vecestim[1:q] <- vecop(Mu0)
vecestim[(q+1):dimens] <- vecop(Sigma0)

# to draw bootstrap samples
#set.seed(2)
bootmatrix <- matrix(sample(n,R*n,replace=TRUE),ncol=R)

bootbiasmat <- matrix(0,dimens,R)  

for (r in 1:R) {
    Yst <- Y[bootmatrix[,r],]
   
    resmatrixst <- Yst - matrix(rep(Mu0,n), n, byrow=TRUE) 
    divecst<- sqrt(mahalanobis(resmatrixst, rep(0,q), Sigma0))
    divecst[divecst<1e-5] <- 1e-5
    udivecst <- rhobiweightder1(divecst,c0)/divecst
    wwdivecst <- rhobiweightder1(divecst,c0)*divecst - rhobiweight(divecst,c0)     
    
    sqrtwu <- sqrt(udivecst)
    B0st <- crossprod(udivecst, Yst) / as.vector(crossprod(sqrtwu))
    wbig <- matrix(rep(sqrtwu,q),ncol=q) 
    wres <- resmatrixst * wbig  
    V0st_term1 <- 1/(b*n)*q * crossprod(wres)
    V0st_term2 <- 1/(b*n)* sum(wwdivecst) * Sigma0
    V0st <- V0st_term1 - V0st_term2
    
    # list uncorrected bootstrap recomputations
    vecfst <- rep(0,dimens)
    vecfst[1:q] <- vecop(B0st)
    vecfst[(q+1):dimens] <- vecop(V0st)
    
    # compute centered, corrected fast bootstrap estimates
    fstbias <- vecfst - vecestim
    bootbiasmat[,r] <- lincorrectmat %*% fstbias  
 
}

##########################################################################

return(list(centered=bootbiasmat, Sest=vecestim))

}
