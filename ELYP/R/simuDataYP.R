simuDataYP <- function(n1=100, n2=100, th1=exp(1), th2=exp(-1), cens=TRUE, alphaX){
#### generate random variables and covariates. Use parameters th1, th2.
#### Two sample case. Sample one is exponential. 
#### Same censoring for both samples. This can improve.
#### alphaX must be of length n1+n2, the first n1 corresponds to Z=0
#### the next n2 corresponds to Z=1. If there is no alpha in the model, then use
#### alphaX= rep(0, n1+n2). Later in the program, we compute exp( alphaX ).
####
#### We use the method of H^(-1)(Z) to generate a r.v. with cumulative
#### hazard H() where Z is a standard exp r.v.

if (length(alphaX) != (n1+n2)) stop("check length of alphaX")

Y1 <- rexp(n1)
EaX1 <- exp( alphaX[1:n1] )
Y1 <- Y1/EaX1

Y2 <- rYP(th1=th1, th2=th2, n=n2, aX=alphaX[(n1+1):(n1+n2)])
Y12 <- c(Y1, Y2)
censor <- rexp(n1+n2, rate=0.5)  ## change for different censor percentage
if(!cens) censor <- censor + Inf
Y <- pmin(Y12, censor)
d <- as.numeric( Y12 <= censor )
Z <- c(rep(0, n1), rep(1,n2))
list(Y=Y, d=d, Zmat=Z)
}
