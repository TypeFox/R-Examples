#' @export
#Exhaustive search over a grid of model parameters
#x: time series
#m: maximum autoregressive order
#d: time delay
#steps: steps ahead

# strategy: instead of running each time a expensive grid search, only first model has grid search, 
# then subsequent use starting values form first model(s) (based on CI as well as point estimates)


selectLSTAR <- function(x, m, d=1, steps=d, mL = 1:m, mH = 1:m, thDelay=0:(m-1), fast=TRUE, trace=FALSE) {
  op <- options(warn=-1)

## matrix of parameters
  IDS <- as.matrix( expand.grid(thDelay, mL, mH) )
  colnames(IDS) <- c("thDelay","mL","mH")

## Loop for each model
  computedAIC <- matrix(NA, nrow=nrow(IDS), ncol=4)
  colnames(computedAIC) <- c("AIC", "BIC", "th", "gamma")
  computedMod <- list()
  for (i in 1:nrow(IDS)){
    mLVal <- IDS[i,2]
    mHVal <- IDS[i,3]
    thDelayVal <- IDS[i,1]
    m <- max(mLVal,mHVal,thDelayVal+1)

  ## for first model compute with full grid
    if(i==1|!fast){# 
      start.list <- list()
  ## for other models, compute with restricted grid
    } else {
      lst <- computedMod[[i-1]]

    ## base values on previous CI
      cf <- confint(lst, parm=c("gamma", "th"), level=0.99) #use previous CI
      if(is.na(cf["gamma",1]) || cf["gamma",1]<0) cf["gamma",1] <-  0.1

    ## base values on previous estimates
      can_th <- sapply(computedMod, getTh) 
      can_gam <- sapply(computedMod, function(x) coef(x)["gamma"])
      can_all <- list(th=can_th, gamma=can_gam)


      start.list <- list(gammaInt=cf["gamma",] , thInt=cf["th",], nGamma=5, nTh=5, candidates=can_all)
    }

    computedMod[[i]]<- lstar(x, m=m, mL=mLVal, mH=mHVal, thDelay=thDelayVal, trace=trace,
		      control=list(maxit=1000), starting.control=start.list)

    if(computedMod[[i]]$model.specific$convergence!=0){
      computedMod[[i]]<- lstar(x, m=m, mL=mLVal, mH=mHVal, thDelay=thDelayVal, trace=trace,
		      control=list(maxit=1000))
    }

  ## for each model, return a few value: AIC, BIC, th and gamma
    computedAIC[i,1] <- AIC( computedMod[[i]])
    computedAIC[i,2] <- BIC( computedMod[[i]])
    computedAIC[i,3] <- getTh( computedMod[[i]])
    computedAIC[i,4] <- coef(computedMod[[i]])["gamma"]

  }


  options(op)

## Sort matrix of parameters according to AIC
  res <- cbind(IDS, computedAIC)
  idSel <- sort(computedAIC[,"AIC"], index=TRUE)$ix
  idSel <- idSel[1:min(10, length(idSel))]
  res <- data.frame(res[idSel,], row.names=NULL)

## Result:
  return(res)
}

if(FALSE){
library(tsDyn)

llynx <- log10(lynx)

environment(lstar) <- environment(star)

system.time(sel1 <- selectLSTAR(x=llynx, m=2))
system.time(sel2_f <- selectLSTAR2(x=llynx, m=2))
system.time(sel2_nof <- selectLSTAR2(x=llynx, m=2, fast=FALSE))


sel1 
sel2_f
sel2_nof


}
