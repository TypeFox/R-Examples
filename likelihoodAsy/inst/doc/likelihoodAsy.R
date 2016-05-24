## ----echo=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(knitr)
opts_chunk$set(comment="", message=FALSE, warning=FALSE,
               tidy.opts=list(keep.blank.line=TRUE, width.cutoff=180),
               fig.width=8, fig.height=6,out.width='1\\textwidth',
               options(width=180),fig.show='hold',fig.align='center',cache=TRUE)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(likelihoodAsy)

## ----Data for Example 1-----------------------------------------------------------------------------------------------------------------------------------------------------------
library(MASS)
data(leuk)

## ----Log likelihood for a Weibull regression model--------------------------------------------------------------------------------------------------------------------------------
loglik.Wbl <- function(theta, data)
{   
  logy <- log(data$y)
  X <- data$X
  loggam <- theta[1]
  beta <- theta[-1]
  gam <- exp(loggam) 
  H <- exp(gam * logy + X %*% beta)
  out <- sum(X %*% beta + loggam + (gam - 1) * logy - H)
  return(out)  
}

## ----Data list for Example 1------------------------------------------------------------------------------------------------------------------------------------------------------
X <- model.matrix(~log(wbc, base=10), data=leuk[leuk$ag=="present",])
data.fz <-list(X = X, y = leuk$time[leuk$ag=="present"])

## ----Data generation for a Weibull regression model-------------------------------------------------------------------------------------------------------------------------------
gendat.Wbl <- function(theta, data)   
{  
  X <- data$X
  n <- nrow(X)
  beta <- theta[-1]
  gam <- exp(theta[1])
  data$y <- (rexp(n) / exp(X %*% beta)) ^ (1 / gam)
  return(data)  
}

## ----Scalar function of interest--------------------------------------------------------------------------------------------------------------------------------------------------
psifcn.Wbl <- function(theta)
{ 
  beta <- theta[-1]
  gam <- exp(theta[1])
  y0 <- 130
  x0 <- 4
  psi <- -(y0 ^ gam) * exp(beta[1] + x0 * beta[2])
  return(psi)   
}

## ----Testing a value for the log Survival function--------------------------------------------------------------------------------------------------------------------------------
rs <- rstar(data=data.fz, thetainit = c(0, 0, 0), floglik = loglik.Wbl, 
            fpsi = psifcn.Wbl, psival = log(0.03), datagen = gendat.Wbl, 
            trace=FALSE, seed=10, psidesc="Log survival function")
rs

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
summary(rs)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
grad.Wbl <- function(theta, data)
{
  logy <- log(data$y)
  X <- data$X
  loggam <- theta[1]
  beta <- theta[-1]
  gam <- exp(loggam)
  H <- exp(gam * logy + X %*% beta)
  score.beta <- t(X) %*% (1 - H) 
  score.nu <- sum(1 + gam * logy - gam * H * logy)
  out <- c(score.nu, score.beta)
  return(out)
}

## ----Checking the gradient--------------------------------------------------------------------------------------------------------------------------------------------------------
cbind(grad(loglik.Wbl, rs$theta.hyp, data=data.fz), 
      grad.Wbl(rs$theta.hyp, data.fz))

## ----Confidence intervals for the log Survival function---------------------------------------------------------------------------------------------------------------------------
rs.int <- rstar.ci(data=data.fz, thetainit = c(0, 0, 0), floglik = loglik.Wbl,
                   fpsi = psifcn.Wbl, fscore=grad.Wbl, datagen=gendat.Wbl, 
                   trace=FALSE, seed=1223, psidesc="Log survival function")
rs.int

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
summary(rs.int)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(rs.int)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print(exp(rs.int$CIrs[2,]), digits=3)

## ----Log likelihood for the AR1 example-------------------------------------------------------------------------------------------------------------------------------------------
likAR1 <- function(theta, data)
{      
  y <- data$y
  mu <- theta[1]
  sigma2 <- exp(theta[2] * 2)
  rho <- theta[3]
  n <- length(y)
  Gamma1 <- diag(1 + c(0, rep(rho^2, n-2), 0))
  for(i in 2:n)
    Gamma1[i,i-1]<- Gamma1[i-1,i] <- -rho 
  lik <- -n/2 * log(sigma2) + 0.5 * log(1 - rho^2) - 1 / (2 * sigma2) * 
          mahalanobis(y, rep(mu,n), Gamma1, inverted = TRUE)
  return(lik)
}

## ----Gradient function for the AR1 example----------------------------------------------------------------------------------------------------------------------------------------
grAR1 <- function(theta, data)
{ 
  y <- data$y
  mu <- theta[1]
  sigma2 <- exp(theta[2] * 2)
  rho <- theta[3]
  n <- length(y)
  Gamma1 <- diag( 1 + c(0, rep(rho^2, n-2), 0))
  DGamma1 <- diag(c(0, rep( 2 * rho, n-2), 0))
  for(i in 2:n)
  {
    Gamma1[i,i-1]<- Gamma1[i-1,i] <- -rho
    DGamma1[i,i-1] <- DGamma1[i-1,i] <- -1   
  }
  out <- rep(0, length(theta))
  out[1] <-  1 / sigma2 * t(rep(1,n)) %*% Gamma1 %*% (y-mu)
  out[2] <-  -n / (2 * sigma2) + 1 / (2 * sigma2^2) * 
             mahalanobis(y, rep(mu,n), Gamma1, inverted = TRUE)
  out[2] <- out[2] * sigma2 * 2
  out[3] <-  -rho / (1 - rho^2) - 1 / (2 * sigma2) * 
             mahalanobis(y, rep(mu,n), DGamma1, inverted = TRUE)
  return(out)
}

## ----Data generation for the AR1 example------------------------------------------------------------------------------------------------------------------------------------------
genDataAR1 <- function(theta, data)  
{
  out <- data
  mu <- theta[1]
  sigma <- exp(theta[2])
  rho <- theta[3]
  n <- length(data$y)
  y <- rep(0,n)
  y[1] <- rnorm(1, mu, s = sigma * sqrt(1 / (1 - rho^2)))
  for(i in 2:n)
    y[i] <- mu + rho * (y[i-1] - mu) + rnorm(1) * sigma 
  out$y <- y
  return(out)
}

## ----Data for the AR1 example-----------------------------------------------------------------------------------------------------------------------------------------------------
data.AR1 <- list( y = as.numeric(lh) )

## ----Testing a value for the correlation------------------------------------------------------------------------------------------------------------------------------------------
rsAR1 <- rstar(data=data.AR1, thetainit = c(0, 0, 0), floglik = likAR1, 
               fpsi = function(theta) theta[3], fscore=grAR1,
            psival = 0.765, datagen=genDataAR1, trace=FALSE, seed=10121,
            psidesc="Autocorrelation parameter")
summary(rsAR1)

## ----Comparison with parametric bootstrap, eval=FALSE-----------------------------------------------------------------------------------------------------------------------------
#  rvals <- rep(0, 10000)
#  set.seed(107)
#  for(i in 1:length(rvals))
#  {
#    data.boot <- genDataAR1(rsAR1$theta.hyp, data.AR1)
#    if(i%%1000==0) cat("i=",i,"\n")
#    r.boot <- rstar(data=data.boot, thetainit = rsAR1$theta.hyp, floglik = likAR1,
#              fpsi = function(theta) theta[3], fscore=grAR1,
#              psival = 0.765, datagen=genDataAR1, trace=FALSE, ronly=TRUE)
#    rvals[i] <- r.boot$r
#  }

## ----Comparison with parametric bootstrap: p-values, eval=FALSE-------------------------------------------------------------------------------------------------------------------
#   c(mean(rvals < rsAR1$r), pnorm(rsAR1$rs) )

## ----Data definition for binomial overdispersion----------------------------------------------------------------------------------------------------------------------------------
data(finndat)
z <- scale(finndat$z * 10, scale=FALSE)
X <- cbind(rep(1,length(z)), z)
data.binOD  <- list(X=X, den = finndat$den, y = finndat$y, 
                    gq=gauss.quad(40,"hermite"))

## ----Log likelihood  function for binomial overdispersion-------------------------------------------------------------------------------------------------------------------------
loglik.binOD <- function(theta, data)
{
  p.range<- function(p, eps=2.22e-15)
   {  
      out <- p
      out[p<eps] <- eps
      out[p>(1-eps)] <- (1-eps)
      return(out)
   }
  y <- data$y 
  den <- data$den    
  X <- data$X  
  gq <- data$gq
  n <- length(y) 
  p <- ncol(X)
  beta <- theta[1:p]  
  sigma <-  exp(theta[p+1])
  linpred <- X %*% beta
  L <- rep(0,n)
  for (i in 1:n)
  {
    prob <- p.range(plogis(linpred[i] + gq$nodes * sqrt(2)*sigma))
    likq <- y[i] * log(prob) + (den[i] - y[i]) * log(1-prob)
    L[i] <- sum(gq$weights * exp(likq) ) / sqrt(2 * pi)
    }
   return(log(prod(L)))
}

## ----Gradient function for binomial overdispersion--------------------------------------------------------------------------------------------------------------------------------
grad.binOD <- function(theta,data)
{ 
  p.range<- function(p, eps=2.22e-15)
   {  
      out <- p
      out[p<eps] <- eps
      out[p>(1-eps)] <- (1-eps)
      return(out)
   }
  y <- data$y 
  den <- data$den    
  X <- data$X  
  gq <- data$gq
  n <- length(y) 
  p <- ncol(X)
  beta <- theta[1:p]  
  sigma <-  exp(theta[p+1])
  linpred <- X %*% beta
  L <- rep(0,n)
  LB <- matrix(0, nrow=n, ncol=p+1)
  out <- rep(0,p+1)
  for (i in 1:n)
  {
    prob <- p.range(plogis(linpred[i]+gq$nodes*sqrt(2)*sigma))
    likq <- y[i] * log(prob) + (den[i] - y[i]) * log(1-prob)
    score <- (y[i] - den[i] * prob)  
    L[i] <- sum(gq$weights * exp(likq) ) / sqrt(2 * pi)
    LB[i,1] <- sum(gq$weights * exp(likq) * score) / sqrt(2 * pi)
    LB[i,2] <- sum(gq$weights * exp(likq) * score * X[i,2] ) / sqrt(2 * pi)
    LB[i,3] <- sum(gq$weights * exp(likq) * score * gq$nodes * 
                  sqrt(2)) / sqrt(2 * pi) * sigma
    out <- out + LB[i,] / L[i]
    }
   return(out)
}

## ----Function that generates a data set for  binomial overdispersion--------------------------------------------------------------------------------------------------------------
gendat.binOD <- function(theta, data)
 {
   out <- data
   den <- data$den
   X <- data$X  
   p <- ncol(X)
   n <- length(data$y)
   beta <- theta[1:p]
   sigma <-  exp(theta[p+1])
   u <- rnorm(n) * sigma
   linpred <- X %*% beta + u
   out$y <- rbinom(n, size=den, prob=plogis(linpred))
   return(out) 
 }

## ----Testing that the slope is one------------------------------------------------------------------------------------------------------------------------------------------------
rs <- rstar(data=data.binOD, thetainit=c(0, 0, 0),  floglik=loglik.binOD, 
            fscore=grad.binOD,  fpsi=function(theta) return(theta[2]), seed=110,
            trace=FALSE, R=500, psival=1 ,datagen=gendat.binOD, 
            psidesc="Regression slope") 
summary(rs)

## ----Log likelihood and data generation for 2x2 table-----------------------------------------------------------------------------------------------------------------------------
loglik.Pois <- function(theta, data)
{
  y <- data$y
  y <- y + 0.50 * c(-1,1,1,-1) ### continuity correction
  mu <- exp(data$X %*% theta)
  el <- sum(y * log(mu) - mu)
  return(el)
}

gendat.Pois <- function(theta, data)
{
  out <- data
  mu <-  exp(data$X %*% theta)
  out$y <- rpois(n=4, lam=mu)
  return(out)
}

## ----Data definition for 2x2 table------------------------------------------------------------------------------------------------------------------------------------------------
rowf <- c(1, 0, 1, 0)
colf <- c(1, 1, 0, 0)
intf <- c(0, 0, 0, 1)
X <- cbind( rep(1, 4), rowf, colf, intf)
data.2x2  <- list(y = c(15, 9, 7, 13), X=X)

## ----Testing independence in 2 x 2 table------------------------------------------------------------------------------------------------------------------------------------------
rs <- rstar(data=data.2x2, thetainit = c(0, 0, 0, 0), floglik = loglik.Pois, 
            fpsi = function(theta) theta[4], psival = 0, datagen=gendat.Pois, 
            trace=FALSE, R=50, psidesc="Independence test")
summary(rs)

## ----Accessing the crying babies data---------------------------------------------------------------------------------------------------------------------------------------------
library(cond)
data(babies)

## ----Standard logistic regression model-------------------------------------------------------------------------------------------------------------------------------------------
mod.glm <- glm(formula = cbind(r1, r2) ~ day + lull - 1, family = binomial, 
               data = babies)
data.obj <- list(y = babies$r1, den = babies$r1 + babies$r2, 
                 X = model.matrix(mod.glm))

## ----Functions for logistic regression--------------------------------------------------------------------------------------------------------------------------------------------
loglik.logit<- function(theta, data) 
{
  y <- data$y
  den <- data$den
  X <- data$X
  eta <- X %*% theta
  p <- plogis(eta)
  l <- sum(y * log(p) + (den - y) * log(1-p))
  return(l)
}

grad.logit<- function(theta, data) 
{
  y <- data$y
  den <- data$den
  X <- data$X
  eta <- X %*% theta
  p <- plogis(eta)
  out <- t(y - p * den) %*% X
  return(drop(out))
}


gendat.logit<- function(theta, data)
{
  X <- data$X
  eta <- X %*% theta
  p <- plogis(eta)
  out <- data
  out$y <- rbinom(length(data$y), size = data$den, prob = p)
  return(out) 
}	

## ----Confidence intervals for the coefficient of lull-----------------------------------------------------------------------------------------------------------------------------
time.with <- system.time( rs.int <- rstar.ci(data=data.obj, 
                         thetainit = coef(mod.glm), 
                         floglik = loglik.logit, fpsi = function(theta) theta[19], 
                         fscore=grad.logit, datagen=gendat.logit, trace=FALSE, 
                          psidesc="Coefficient of lull") ) 
time.without <- system.time( rs.int.no <- rstar.ci(data=data.obj, 
                         thetainit = coef(mod.glm), 
                         floglik = loglik.logit, fpsi = function(theta) theta[19], 
                         datagen=gendat.logit, trace=FALSE, 
                         psidesc="Coefficient of lull") )

## ----Summary of confidence intervals for the coefficient of lull------------------------------------------------------------------------------------------------------------------
summary(rs.int)

## ----Comparison of computational times--------------------------------------------------------------------------------------------------------------------------------------------
time.with
time.without

## ----Results obtained with the cond package---------------------------------------------------------------------------------------------------------------------------------------
res.cond <- cond(object = mod.glm, offset = lullyes)
summary(res.cond)

## ----Numerical optimization of profile and modified profile log likelihoods-------------------------------------------------------------------------------------------------------
max.prof <- nlminb(0, logPL, data=data.obj, thetainit=coef(mod.glm), 
                  floglik=loglik.logit, fscore=grad.logit, indpsi=19, trace=FALSE, 
                  minus=TRUE)
max.mpl <- nlminb(0, logMPL, data=data.obj, mle=coef(mod.glm), 
                  floglik=loglik.logit, fscore=grad.logit, datagen=gendat.logit,
                  indpsi=19, R=50, seed=2020, trace=FALSE, minus=TRUE)
c(max.prof$par, max.mpl$par)

## ----Plotting the two log likelhoods----------------------------------------------------------------------------------------------------------------------------------------------
psi.vals <- seq(-0.3, 3.7, l=30)
obj.prof <- sapply(psi.vals, logPL, data=data.obj, thetainit=coef(mod.glm), 
                floglik=loglik.logit, fscore=grad.logit, indpsi=19)
obj.mpl <- sapply(psi.vals, logMPL, data=data.obj, mle=coef(mod.glm), 
                floglik=loglik.logit, fscore=grad.logit, datagen=gendat.logit,
                indpsi=19, R=50, seed=2020)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
par(pch="s")
plot(psi.vals, obj.prof - max(obj.prof), type="l", xlab=expression(psi), 
     ylab="log likelihood", lwd=2, las=1)
lines(psi.vals, obj.mpl  - max(obj.mpl), col="red", lwd=2)
legend("topright", col=c(1, 2), lty=1, lwd=2, legend=c("Profile","MPL"), bty="n")

## ----Functions for random intercept model-----------------------------------------------------------------------------------------------------------------------------------------
logLikLme<- function(theta, data) 
{   
  X <- data$X
  Z <- data$Z
  y <- data$y
  beta <- theta[1:ncol(X)]
  sigma.b <- theta[ncol(X)+1]
  sigma <- theta[ncol(X)+2]
  n <- nrow(X)
  V <- tcrossprod(Z) * sigma.b^2 + diag(n) * sigma^2 
  L <- chol(V)
  XL <- backsolve(L, X, transpose=TRUE)
  yL <- backsolve(L, y, transpose=TRUE)
  out<- - sum(log(diag(L))) - sum( (yL-XL %*% beta)^2) / 2
  return(out)
}


gradLikLme <- function(theta, data) 
{   
  X <- data$X
  Z <- data$Z
  y <- data$y
  beta <- theta[1:ncol(X)]
  sigma.b <- theta[ncol(X)+1]
  sigma <- theta[ncol(X)+2]
  n <- nrow(X)
  V <- tcrossprod(Z) * sigma.b^2 + diag(n) * sigma^2     
  L <- chol(V)
  XL<- backsolve(L, X, transpose=TRUE)   
  yL<- backsolve(L, y, transpose=TRUE)
  out <- rep(0, length(theta))
  out[1:ncol(X)] <-  t(yL-XL %*% beta)  %*% XL
  ni<- as.vector(t(Z) %*% rep(1,n))
  Zv<- matvec(Z, sqrt(1/(sigma^2 + sigma.b^2 * ni)))
  V1 <- diag(n) / sigma^2 - tcrossprod(Zv) * sigma.b^2 / sigma^2
  Vb <- tcrossprod(Z) * 2 * sigma.b
  Vs <- diag(n) * 2 * sigma
  Mb <- V1 %*% Vb
  Ms <- V1 %*% Vs
  r <- as.vector(y - X %*% beta)
  out[ncol(X)+1] <- -sum(diag(Mb)) / 2  + 
                    as.numeric( t(r) %*% Mb %*% V1 %*% r) / 2  
  out[ncol(X)+2] <- -sum(diag(Ms)) / 2  + 
                    as.numeric( t(r) %*% Ms %*% V1 %*% r) / 2 
  return(out)
}


genDataLme <- function(theta, data)   
{
  out <- data
  X <- data$X
  Z <- data$Z
  y <- data$y
  beta <- theta[1:ncol(X)]
  sigma.b <- theta[ncol(X)+1]
  sigma <- theta[ncol(X)+2]
  n <- nrow(X)
  mu <- X %*% beta
  b <- rnorm(ncol(Z), s=sigma.b)
  e <- rnorm(nrow(Z), s=sigma)
  out$y <- mu + e + Z %*% b
  return(out)
}

## ----Random intercept example-----------------------------------------------------------------------------------------------------------------------------------------------------
library(lme4)
fm1R <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
sleepdata <- list(X=model.matrix(Reaction ~ Days, sleepstudy),
                  Z=model.matrix(Reaction ~ factor(Subject)-1, sleepstudy),
                  y=sleepstudy$Reaction)

## ----MLE for example--------------------------------------------------------------------------------------------------------------------------------------------------------------
mleFull <- optim( c(250, 10, 30, 30), logLikLme, gr=gradLikLme, 
                  data=sleepdata, method="BFGS",
                  control=list(fnscale=-1)) 

## ----Modified Profile likelihood maximization, eval=FALSE-------------------------------------------------------------------------------------------------------------------------
#  mleM <- optim(mleFull$par[3:4], logMPL, data=sleepdata, mle=mleFull$par,
#                floglik=logLikLme, fscore=gradLikLme, minus=TRUE,
#                indpsi=3:4, datagen=genDataLme, trace=FALSE, seed=11, R=100)

