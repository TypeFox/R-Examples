Sboot_twosample <-function(X, groups, R=999, ests=Sest_twosample(X,groups))
{
# Robust bootstrap for two populations S-estimator
# INPUT:
# X = data matrix
# R = number of bootstrap samples
# ests = result of twosampleS
# OUTPUT:
# centered = bootstrap recalculations (centered by original estimates)
#            = (p+p+p*p) x R  (all parameters are stacked into one column)
# Sest = ((p+p + p*p) x 1) original two population S estimates in vec-form
#

# 

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



#------------------------------------------------------------------------ 
#-                                main function                         -
#------------------------------------------------------------------------
n <- nrow(X)
p <- ncol(X)

n1 <- sum(groups==1)
n2 <- sum(groups==2)
X1 <- X[groups==1,,drop=FALSE]     
X2 <- X[groups==2,,drop=FALSE]
dimens <- 2*p+p*p

Ip <- diag(p)

# determine the constants in Tukey biweight S:
c <- ests$c
b <- ests$b

mu1 <- ests$Mu1
mu2 <- ests$Mu2
Sigma <- ests$Sigma
Sinv <- solve(Sigma)

###############################################################################
###                calculate jacobian for correction matrix                 ###
###############################################################################

Resi1 <- X1 - matrix(rep(mu1,n1), n1, byrow=TRUE)
divec1 <- sqrt(mahalanobis(Resi1, rep(0,p), Sigma))
divec1[divec1<1e-5] <- 1e-5

udivec1 <- rhobiweightder1(divec1,c)/divec1
wdivec1 <- (rhobiweightder2(divec1,c)*divec1 - rhobiweightder1(divec1,c))/ divec1^3
zdivec1 <- rhobiweightder2(divec1,c)
wwdivec1 <- rhobiweightder1(divec1,c)*divec1-rhobiweight(divec1,c)

Resi2 <- X2 - matrix(rep(mu2,n2), n2, byrow=TRUE) 
divec2 <- sqrt(mahalanobis(Resi2, rep(0,p), Sigma))
divec2[divec2<1e-5] <- 1e-5

udivec2 <- rhobiweightder1(divec2,c)/divec2
wdivec2 <- (rhobiweightder2(divec2,c)*divec2 - rhobiweightder1(divec2,c))/divec2^3
zdivec2 <- rhobiweightder2(divec2,c)
wwdivec2 <- rhobiweightder1(divec2,c)*divec2 - rhobiweight(divec2,c)

an1 <- sum(udivec1)
Vn1 <- crossprod(udivec1,X1)
an2 <- sum(udivec2)
Vn2 <- crossprod(udivec2,X2)

term1a <- matrix(0,1,p) 
term1b <- matrix(0,p,p)
term3a <- matrix(0,1,p*p) 
term3b <- matrix(0,p,p*p)
term5a <- matrix(0,1,p) 
term5b <- matrix(0,p,p)
term6a <- matrix(0,1,p*p) 
term6b <- matrix(0,p,p*p)
term7a <- matrix(0,p*p,p) 
term7b <- matrix(0,p*p,p) 
term7c <- matrix(0,1,p)
term8a <- matrix(0,p*p,p) 
term8b <- matrix(0,p*p,p) 
term8c <- matrix(0,1,p)
term9a <- matrix(0,p*p,p*p) 
term9b <- matrix(0,1,p*p) 
term9c <- matrix(0,p*p,p*p) 
term9d <- matrix(0,1,p*p)

for (i in 1:n1) {
    Xi <- as.matrix(X1[i,])
    resi <- as.matrix(Resi1[i,])
    vecresiresi <- vecop(resi%*%t(resi))
    wdi <- wdivec1[i]
    zdi <- zdivec1[i]
    udi <- udivec1[i]
    tveci <- crossprod(resi,Sinv)
    tvecSi <- t(vecop(Sinv%*%resi%*%t(resi)%*%Sinv))
    
    term1a <- term1a + wdi*t(resi)
    term1b <- term1b + wdi*(Xi%*%tveci)

    term3a <- term3a + wdi*tvecSi
    term3b <- term3b + wdi*(Xi%*%tvecSi)

    term7a <- term7a + wdi*(vecresiresi%*%tveci)
    term7b <- term7b + udi*(kronecker(Ip,resi)+kronecker(resi,Ip))
    term7c <- term7c + zdi*tveci
    
    term9a <- term9a + wdi*(vecresiresi%*%tvecSi)
    term9b <- term9b + zdi*tvecSi
}
for (i in 1:n2) {
    Xi <- as.matrix(X2[i,])
    resi <- as.matrix(Resi2[i,])
    vecresiresi <- vecop(resi%*%t(resi))
    wdi <- wdivec2[i]
    zdi <- zdivec2[i]
    udi <- udivec2[i]
    tveci <- crossprod(resi,Sinv)
    tvecSi <- t(vecop(Sinv%*%resi%*%t(resi)%*%Sinv))
    
    term5a <- term5a + wdi*t(resi)
    term5b <- term5b + wdi*(Xi%*%tveci)

    term6a <- term6a + wdi*tvecSi
    term6b <- term6b + wdi*(Xi%*%tvecSi)

    term8a <- term8a + wdi*(vecresiresi%*%tveci)
    term8b <- term8b + udi*(kronecker(Ip,resi)+kronecker(resi,Ip))
    term8c <- term8c + zdi*tveci
    
    term9c <- term9c + wdi*(vecresiresi%*%tvecSi)
    term9d <- term9d + zdi*tvecSi
}

partder1 <- as.matrix(t(Vn1)/an1^2)%*%term1a%*%Sinv - term1b/an1
partder2 <- matrix(0,p,p)
partder3 <- 1/2*as.matrix(t(Vn1)/an1^2)%*%term3a - 1/2*term3b/an1
partder4 <- matrix(0,p,p)
partder5 <- as.matrix(t(Vn2)/an2^2)%*%term5a%*%Sinv - term5b/an2
partder6 <- 1/2*as.matrix(t(Vn2)/an2^2)%*%term6a - 1/2*term6b/an2
partder7 <- -p/(b*n)*(term7a + term7b) + 1/(b*n)*vecop(Sigma)%*%term7c
partder8 <- -p/(b*n)*(term8a + term8b) + 1/(b*n)*vecop(Sigma)%*%term8c
partder9 <- -p/(2*b*n)*(term9a + term9c) + 1/(2*b*n)*vecop(Sigma)%*%(term9b+term9d) - 1/(b*n)*(sum(wwdivec1)+sum(wwdivec2))*diag(rep(1,p*p))


jacobian <- diag(rep(0,dimens))
jacobian[1:p,1:p] <- partder1
jacobian[1:p,(p+1):(2*p)] <- partder2
jacobian[1:p,(2*p+1):dimens] <- partder3
jacobian[(p+1):(2*p),1:p] <- partder4
jacobian[(p+1):(2*p),(p+1):(2*p)] <- partder5
jacobian[(p+1):(2*p),(2*p+1):dimens] <- partder6
jacobian[(2*p+1):dimens,1:p] <- partder7
jacobian[(2*p+1):dimens,(p+1):(2*p)] <- partder8
jacobian[(2*p+1):dimens,(2*p+1):dimens] <- partder9

Idim <- diag(dimens)
lincorrectmat <- solve(Idim-jacobian)

##############################################################################

# now do the actual bootstrapping

# put all parameters into one large vector:
vecestim <- rep(0,dimens)
vecestim[1:p] <- mu1
vecestim[(p+1):(2*p)] <- mu2
vecestim[(2*p+1):dimens] <- vecop(Sigma)

# draw bootstrap samples for each group separately:
set.seed(2)
bootmatrix1 <- matrix(sample(n1,R*n1,replace=TRUE),ncol=R)
bootmatrix2 <- matrix(sample(n2,R*n2,replace=TRUE),ncol=R)

bootbiasmat <- matrix(0,dimens,R)  
bootbiaszc <- matrix(0,dimens,R) 

for (r in 1:R) {
    X1st <- X1[bootmatrix1[,r],,drop=FALSE]
    Resi1st <- X1st - matrix(rep(mu1,n1), n1, byrow=TRUE)
    divec1st <- sqrt(mahalanobis(Resi1st, rep(0,p), Sigma))
    divec1st[divec1st<1e-5] <- 1e-5
    
    udivec1st <- rhobiweightder1(divec1st,c)/divec1st

    wwdivec1st <- rhobiweightder1(divec1st,c)*divec1st - rhobiweight(divec1st,c)   

    X2st <- X2[bootmatrix2[,r],,drop=FALSE]
    Resi2st <- X2st -matrix(rep(mu2,n2), n2, byrow=TRUE) 
    divec2st <- sqrt(mahalanobis(Resi2st, rep(0,p), Sigma))
    divec2st[divec2st<1e-5] <- 1e-5
    
    udivec2st <- rhobiweightder1(divec2st,c)/divec2st
    wwdivec2st <- rhobiweightder1(divec2st,c)*divec2st - rhobiweight(divec2st,c)   
    
    an1st <- sum(udivec1st)
    Vn1st <- crossprod(udivec1st,X1st)
    sqrtwu1 <- sqrt(udivec1st)
    wbig1 <- matrix(rep(sqrtwu1,p),ncol=p) 
    wres1 <- Resi1st * wbig1  
        
    Bn1st <- p* crossprod(wres1)
    
    an2st <- sum(udivec2st)
    Vn2st <- crossprod(udivec2st,X2st)
    sqrtwu2 <- sqrt(udivec2st)
    wbig2 <- matrix(rep(sqrtwu2,p),ncol=p) 
    wres2 <- Resi2st * wbig2  
        
    Bn2st <- p* crossprod(wres2)
    
    
    vecfst <- rep(0,dimens)
    vecfst[1:p] <- (1/an1st)*Vn1st
    vecfst[(p+1):(2*p)] <- (1/an2st)*Vn2st

    Sigmast <- 1/(b*n)*(Bn1st + Bn2st) - 1/(b*n)*(sum(wwdivec1st)+ sum(wwdivec2st)) * Sigma
    vecfst[(2*p+1):dimens] <- vecop(Sigmast)
    
	
    fstbias <- vecfst - vecestim
    
    bootbiaszc[,r] <- fstbias
    
    # apply linear correction:
    bootbiasmat[,r] <- lincorrectmat %*% fstbias  
}    

return(list(centered=bootbiasmat,Sest=vecestim))
}

###############################################################################

