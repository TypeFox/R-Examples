rwz.reject <-
function(Y=NULL, nlags=4, draws=200, subdraws=200, nkeep=1000, KMIN=1, KMAX=4, constrained=NULL, constant=TRUE, steps=24){
#
#--- SANITY CHECK ---#
sanity.check.reject(Y=Y, nlags=nlags, draws=draws, subdraws=subdraws, nkeep=nkeep, KMIN=KMIN, KMAX=KMAX, constrained=constrained, constant=constant, steps=steps)
#
#--- SET UP PARAS ---#
varnames <- colnames(Y)
n1 <- draws
n2 <- subdraws
nstep <- steps
nlags <- nlags
nvar <- ncol(Y)
nobs <- nrow(Y)
nnobs0 <- nlags + 1
nnobs <- nobs - nlags
nnvar0 <- nvar + 1
ntot <- n1*n2
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
#
#--- SIGMA and SXX ---#
uu <- crossprod(resid)
# sigma <- (1/(nnobs-nncoef))*uu
sigma <- (1/nnobs)*uu
#
#--- SET UP MCMC OF VAR ---#
sxx <-  chol(xx)
sv <- solve(uu)
svt <-  chol(sv)
betaols <- t(bcoef)
best <- betaols
wishdof <- nnobs-nncoef
#
#--- MATRICES FOR DRAWS ---#
goodresp <- array(NA, c(nkeep, nstep, nvar))
BDraws <- array(NA, c(n1, nncoef, nvar))
SDraws <- array(NA, c(n1, nvar, nvar))
imp <- matrix(NA, nrow=nstep, ncol=nvar)
fevd <- matrix(NA, nrow=nstep, ncol=nvar)
goodfevd <- array(NA, c(nkeep, nstep, nvar))
goodshock <- array(NA, c(nkeep, nnobs))
uhatt <- matrix(NA, nnobs, 1)
#
#--- MCMC INTEGRATION ---#
accept <- 0
message('Starting MCMC, ', date(),'.', sep="")
pb0 <- txtProgressBar(min = 0, max = n1, style = 3)
for(draws in 1:n1){
  setTxtProgressBar(pb0, draws)
  #
  #--- sigma draws ---#
  sigmad  <- solve(matrix(rWishart(1, wishdof, sv), nrow=nvar, ncol=nvar))
  swish   <- chol(sigmad)
  #
  #--- beta draws ---#
  swsxx <-   sigmad  %x% xx
  bd <- rep(0, nrow(swsxx))
  # betau <- matrix(mvrnormR(1,0,swsxx), nrow=nncoef, ncol=nvar)
   betau <- matrix(mvnfast::rmvn(1, bd, swsxx), nrow=nncoef, ncol=nvar)
  betadraw <- betaols + betau
  bhat <- betadraw
  #
  #--- irfs ---#
  imfhat <- fn.impulse(bhat, swish, c(nvar, nlags, nstep))
  impulses <-  array(imfhat, dim=c(nstep,nvar,nvar))
  imp2 <- impulses^2
  imp2sum <- apply(imp2, c(2,3), cumsum)
  mse <-  apply(imp2sum, c(1,2), sum)
  fevd0 <- array(apply(imp2sum, 3, "/",  mse), dim=c(nstep, nvar, nvar))
  #
  for(subdraws in 1:n2){
    a <- matrix(rnorm(nvar,mean=0,sd=1), nvar, 1)
    RWZA <- RWZAccept(a,KMIN,KMAX, constrained, impulses)
    RWZ <- RWZA$acc
    q <- RWZA$Q
    #
    if(RWZ==1){
      for(j in 1:nstep){
        imp[j,] <- t(impulses[j,,]%*%q)
        fevd[j,] <- t(fevd0[j,,]%*%(q^2))
      }
      accept <- accept+1
      goodresp[accept, ,] <-  imp
      goodfevd[accept, ,] <- fevd * 100
      BDraws[draws, , ] <- betadraw
      SDraws[draws, , ] <- sigmad
      uhat <-   Y[nnobs0:nobs ,] - data %*% bhat
      for(i in 1:nnobs){
        uhatt[i,] <-   uhat[i, ] %*%  (  solve(swish) %*% q)
      }
      goodshock[accept, ] <-  t(uhatt)
    }else{
      next
    }
    if(accept>=nkeep){
      break
    }
  } # end subdraws
#
if(accept>=nkeep){
  break
}
  ldraw <- draws
}#END DRAWS
close(pb0)
#
#--- FIX PARA MATRICES ---#
if(ldraw<n1){
  BDraws <- BDraws[1:ldraw, , ]
  SDraws <- SDraws[1:ldraw, , ]
  dimnames(SDraws) <- list(1:ldraw, varnames, varnames)
}
#
#--- WARNING MESSAGE IN CASE OF TOO FEW DRAWS ---#
if(accept<nkeep){#this can be shortened
  if(accept==0){
    stop("\n Not enough accepted draws to proceed!")
  }else{
    goodresp <- goodresp[1:accept, , ]
    goodfevd <- goodfevd[1:accept, , ]
    goodshock <- goodshock[1:accept, ]
    message('\n Warning! Had only ', accept,' accepted draw(s) out of ',ntot,'.', sep="")
  }
}
nn1 <- accept
dimnames(goodresp) <- list(1:nn1, 1:nstep, varnames)
dimnames(goodfevd) <- list(1:nn1, 1:nstep, varnames)
#
if(constant == FALSE){
  dimnames(BDraws) <-  list(1:ldraw, c(paste(varnames,rep(1:nlags, each=length(varnames)), sep="")) , varnames)}else{
    dimnames(BDraws) <- list(1:ldraw, c(paste(varnames,rep(1:nlags, each=length(varnames)), sep=""),"const"), varnames)
  }
#
message('\n MCMC finished, ', date(),'.', sep="")
return(list(IRFS=goodresp, FEVDS = goodfevd,  SHOCKS = goodshock, BDraws=BDraws, SDraws=SDraws))
}
