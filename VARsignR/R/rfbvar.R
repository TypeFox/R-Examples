rfbvar <-
function(Y=NULL, nlags=4, draws=1000, constant=TRUE, steps=24, shock=1){
#
#--- SANITY CHECK ---#
sanity.check.rfbvar(Y=Y, nlags=nlags, draws=draws, constant=constant, steps=steps, shock=shock)

#--- SET UP PARAS ---#
varnames <- colnames(Y)
n1 <- draws
nstep <- steps
nlags <- nlags
shock <- shock
nvar <- ncol(Y)
nobs <- nrow(Y)
nnobs0 <- nlags + 1
nnobs <- nobs - nlags
nnvar0 <- nvar + 1
#
if(constant == FALSE){
CONS <- "F"
ncoef <- nvar * nlags
nncoef <- nvar * nlags
nnvar1 <- nvar * (nlags + 1)
}else{
CONS <- "T"
ncoef <- nvar * (nlags+1)
nncoef <- nvar * nlags + 1
nnvar1 <- nvar * (nlags + 1) + 1
}
#
#---REDUCED FORM VAR MODEL ---#
model <- rfvar(Y,lags=nlags, const=CONS)
bcoef <- model$By # same order as above but w/const and nvar x nvar x lags
resid <- model$u # same as above
data <- model$X
xx <- model$xx

#--- SIGMA and SXX ---#
uu <- crossprod(resid)
# sigma <- (1/(nnobs-nncoef))*uu
sigma <- (1/nnobs)*uu

#--- SET UP MCMC OF VAR ---#
sxx <-  chol(xx)
sv <- solve(uu)
svt <-  chol(sv)
betaols <- t(bcoef)
best <- betaols
wishdof <- nnobs-nncoef

#--- MATRICES FOR DRAWS ---#
goodresp <- array(NA, c(n1, nstep, nvar))
BDraws <- array(NA, c(n1, nncoef, nvar))
SDraws <- array(NA, c(n1, nvar, nvar))
imp <- matrix(NA, nrow=nstep, ncol=nvar)
fevd <- matrix(NA, nrow=nstep, ncol=nvar)
goodfevd <- array(NA, c(n1, nstep, nvar))
goodshock <- array(NA, c(n1, nnobs))

#--- Monte CARLO INTEGRATION ---#
message('\n Starting MCMC, ', date(),'.', sep="")
pb <- txtProgressBar(min = 0, max = n1, style = 3)
for(draws in 1:n1){
  setTxtProgressBar(pb, draws)

  #--- sigma draws ---#
  sigmad  <- solve(matrix(rWishart(1, wishdof, sv), nrow=nvar, ncol=nvar))
  swish   <- chol(sigmad)

  #--- beta draws ---#
  swsxx <- sigmad  %x% xx
  bd <- rep(0, nrow(swsxx))
  #betau <- matrix(mvrnormR(1,0,swsxx), nrow=nncoef, ncol=nvar)
  betau <- matrix(mvnfast::rmvn(1, bd, swsxx), nrow=nncoef, ncol=nvar)
  betadraw <- betaols + betau
  bhat <- betadraw

  #--- irfs ---#
  imfhat <- fn.impulse(bhat, swish, c(nvar, nlags, nstep))
  impulses <-  array(imfhat, dim=c(nstep,nvar,nvar))
  imp2 <- impulses^2

  imp2sum <- apply(imp2, c(2,3), cumsum)
  mse <-  apply(imp2sum, c(1,2), sum)
  fevd0 <- array(apply(imp2sum, 3, "/",  mse), dim=c(nstep, nvar, nvar))

  imp <- impulses[,,shock]
  fevd <- fevd0[,,shock]

  goodresp[draws, ,] <- imp
  goodfevd[draws, ,] <- fevd * 100
  BDraws[draws, , ] <- betadraw
  SDraws[draws, , ] <- sigmad
  uhat <-   Y[nnobs0:nobs ,] - data %*%bhat
   uhatt <- t(uhat %*% solve(swish))
   uhatt <-  uhatt[ abs(shock),]
  goodshock[draws, ] <-  uhatt
}#end draws
close(pb)
#
dimnames(goodresp) <- list(1:n1, 1:nstep, varnames)
dimnames(goodfevd) <- list(1:n1, 1:nstep, varnames)
dimnames(SDraws) <- list(1:n1, varnames, varnames)
#
if(constant == FALSE){
dimnames(BDraws) <-  list(1:n1, c(paste(varnames,rep(1:nlags, each=length(varnames)), sep="")) , varnames)}else{
dimnames(BDraws) <- list(1:n1, c(paste(varnames,rep(1:nlags, each=length(varnames)), sep=""),"const"), varnames)
}
#
message('\n MCMC finished, ', date(),'.', sep="")
return(list(IRFS=goodresp, FEVDS =goodfevd, BDraws=BDraws, SDraws=SDraws, SHOCKS = goodshock))
}
