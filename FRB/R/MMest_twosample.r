MMest_twosample<-function(X, groups, control=MMcontrol(...), ...)
{
# M-step for two-sample location and common scatter with auxiliary S-scale
# INPUT:
# X = data matrix
# groups = vector of 1's and 2's, indicating group numbers
# Snsamp = number of random subsamples in fast-S algorithm
# shapeEff = logical; if TRUE, estimate is tuned to 95% shape-efficiency, 
#                                           otherwise 95% location efficiency

# OUTPUT:
# result$Mu1 = MM-estimate first center
# result$Mu2 = MM-estimate second center
# result$Gamma = MM-estimate shape matrix
# result$Sigma = MM-estimate covariance matrix
# result$SMu1 = S-estimate first center
# result$SMu2 = S-estimate second center
# result$SGamma = S-estimate shape matrix
# result$scale = S-estimate scale
# result$SSigma = S-estimate covariance matrix

# --------------------------------------------------------------------

rhobiweight <- function(x,c)
{
# Computes Tukey's biweight rho function with constant c for all values in x

hulp <- x^2/2 - x^4/(2*c^2) + x^6/(6*c^4)
rho <- hulp*(abs(x)<c) + c^2/6*(abs(x)>=c)

return(rho)
}

# --------------------------------------------------------------------

scaledpsibiweight <- function(x,c)
{
# Computes scaled Tukey's biweight psi function with constant c for all values in x

hulp <- 1 - 2*x^2/(c^2) + x^4/(c^4)
psi <- hulp*(abs(x)<c)

return(psi)
}

# --------------------------------------------------------------------
# (taken from Claudia Becker's Sstart0 program)

"chi.int" <- function(p, a, c1)
return(exp(lgamma((p + a)       
        #   partial expectation d in (0,c1) of d^a under chi-squared p
  /2) - lgamma(p/2)) * 2^{a/2} * pchisq(c1^2, p + a))

# --------------------------------------------------------------------

"sigma1.bw" <- function(p, c1)
{  
# called by csolve.bw.MM()

gamma1.1 <- chi.int(p,2,c1) - 6*chi.int(p,4,c1)/(c1^2) + 5*chi.int(p,6,c1)/(c1^4)   
gamma1.2 <- chi.int(p,2,c1) - 2*chi.int(p,4,c1)/(c1^2) + chi.int(p,6,c1)/(c1^4)
gamma1 <- ( gamma1.1 + (p+1)*gamma1.2 ) / (p+2)

sigma1.0 <- chi.int(p,4,c1) - 4*chi.int(p,6,c1)/(c1^2) + 6*chi.int(p,8,c1)/(c1^4) - 4*chi.int(p,10,c1)/(c1^6) + chi.int(p,12,c1)/(c1^8)
return( sigma1.0 / (gamma1^2) * p/(p+2) )

}

# --------------------------------------------------------------------

"loceff.bw" <- function(p, c1)
{  
# called by csolve.bw.MM(); computes location efficiency corresponding to c1

alpha1 <- 1/p * (chi.int(p,2,c1) - 4*chi.int(p,4,c1)/(c1^2) + 6*chi.int(p,6,c1)/(c1^4) - 4*chi.int(p,8,c1)/(c1^6) + chi.int(p,10,c1)/(c1^8))   
beta1.1 <- chi.int(p,0,c1) - 2*chi.int(p,2,c1)/(c1^2) + chi.int(p,4,c1)/(c1^4)
beta1.2 <- chi.int(p,0,c1) - 6*chi.int(p,2,c1)/(c1^2) + 5*chi.int(p,4,c1)/(c1^4)
beta1 <- (1-1/p)*beta1.1 + 1/p*beta1.2

return( beta1^2 / alpha1 )

}

# --------------------------------------------------------------------

csolve.bw.MM <- function(p, eff, shape = TRUE)
{
# constant for second Tukey Biweight rho-function for MM, for fixed shape-efficiency 

maxit <- 1000
eps <- 10^(-8)
diff <- 10^6
#ctest <- csolve.bw.asymp(p,.5)
ctest <- -.4024 + 2.2539 * sqrt(p) # very precise approximation for c corresponding to 50% bdp
iter <- 1
while ((diff>eps) & (iter<maxit)) {
    cold <- ctest
    if (shape == TRUE)
        ctest <- cold * eff * sigma1.bw(p,cold)
    else
        ctest <- cold * eff / loceff.bw(p,cold)
      
    diff <- abs(cold-ctest)
    iter <- iter+1
}
return(ctest)

}


#-------------------------------------------------------------------------
#-                          main function                                -
#-------------------------------------------------------------------------
X <- as.matrix(X)
n <- nrow(X)
p <- ncol(X)

X1 <- X[groups==1,,drop=FALSE]
X2 <- X[groups==2,,drop=FALSE]
n1 <- sum(groups==1)
n2 <- sum(groups==2)

eff <- control$eff
bdp <- control$bdp
shapeEff <- control$shapeEff
fastScontrols <- control$fastScontrols
maxiter <- control$maxIt.MM
mtol <- control$convTol.MM

c <- csolve.bw.MM(p, eff, shape=shapeEff)

# first compute 50% breakdown S-estimator
Sresult <- Sest_twosample(X, groups, bdp, fastScontrols)

c0 <- Sresult$c
b <- Sresult$b

auxscale <- Sresult$scale
newG <- Sresult$Gamma
newmu1 <- Sresult$Mu1
newmu2 <- Sresult$Mu2

R1 <- X1-matrix(rep(newmu1,n1), n1, byrow=TRUE)
R2 <- X2-matrix(rep(newmu2,n2), n2, byrow=TRUE)
Res = rbind(R1,R2)
psres <- sqrt(mahalanobis(Res,rep(0,p),newG))
newobj <- mean(rhobiweight(psres/auxscale,c))
origobj <- newobj

# compute M-estimate with auxilliary scale through IRLS steps, starting
# from S-estimate
iteration <- 1
oldobj <- newobj + 1
while (((oldobj - newobj) > mtol) & (iteration < maxiter)) {
    oldobj <- newobj
    w <- scaledpsibiweight(psres/auxscale,c)
    sqrtw <- sqrt(w)
    mu1new <- crossprod(w[1:n1], X1) / as.vector(crossprod(sqrtw[1:n1]))
    mu2new <- crossprod(w[(n1+1):(n1+n2)], X2) / as.vector(crossprod(sqrtw[(n1+1):(n1+n2)]))
    wbig <- matrix(rep(sqrtw,p),ncol=p) 
    wRes <- Res * wbig  
    newGamma <- crossprod(wRes)
    newGamma <- det(newGamma)^(-1/p)*newGamma
    
    R1new <- X1-matrix(rep(mu1new,n1), n1,byrow=TRUE)
    R2new <- X2-matrix(rep(mu2new,n2), n2,byrow=TRUE)
    Res <- rbind(R1new,R2new)
        
    psres <- sqrt(mahalanobis(Res,rep(0,p),newGamma))
    
    newobj <- mean(rhobiweight(psres/auxscale,c))
    iteration <- iteration+1
}

if (newobj <= origobj) {
    resloc1 <- mu1new
    resloc2 <- mu2new
    resshape <- newGamma
    rescovariance <- newGamma*auxscale^2}
else # isn't supposed to happen
   {resloc1 <- Sresult$loc1
   resloc2 <- Sresult$loc2
   resshape <- Sresult$Gamma
   rescovariance <- Sresult$cov
   }

R1 <- X1 - matrix(rep(resloc1,n1), n1, byrow=TRUE)
R2 <- X2 - matrix(rep(resloc2,n2), n2, byrow=TRUE)
psres <- sqrt(mahalanobis(rbind(R1,R2),rep(0,p),rescovariance))
w <- scaledpsibiweight(psres,c)
outFlag <- (psres > sqrt(qchisq(.975, p)))

return(list(Mu1=resloc1,Mu2=resloc2,Gamma=resshape,Sigma=rescovariance,SMu1=Sresult$Mu1,SMu2=Sresult$Mu2,SGamma=Sresult$Gamma,scale=auxscale,SSigma=Sresult$Sigma,c0=c0,b=b,c1=c,w=w,outFlag=outFlag))
}

