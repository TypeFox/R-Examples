library(survival)

persistence.prob <- function(turbineID=NA, perstime, status, pers.const=FALSE, R=10000){
# dat contains the variables
#            turbineID: name of the turbines
#            perstime: observed persistence times
#            status:  1 = removal was observed, 0 = removal was not observed (the object has still been there at the end of the observations)
# pers.const: set to TRUE if you want to assume constant persistence probability over time. If TRUE an exponential model is used, if FALSE the cox proportional hazard model is used.
# R:          the number or Monte Carlo simulations used to obtain confidence intervals for the estimated persistence probability in the exponential model
#---------------------------------------------------------------------------------
dat <- data.frame(turbineID=factor(turbineID), perstime=perstime, status=status)
newdat <- data.frame(turbineID=factor(levels(dat$turbineID), levels=levels(dat$turbineID)))
if(!pers.const){
  if(nlevels(dat$turbineID)>1){
    surv.mod<-coxph(Surv(perstime, status)~turbineID, data=dat)
    predvalues <- summary(survfit(surv.mod, newdata=newdat))
    estpers <- predvalues$surv
    colnames(estpers) <- levels(dat$turbineID)
    rownames(estpers) <- predvalues$time
    estpers.lwr <- predvalues$lower
    colnames(estpers.lwr) <- levels(dat$turbineID)
    rownames(estpers.lwr) <- predvalues$time
    estpers.upr <- predvalues$upper
    colnames(estpers.upr) <- levels(dat$turbineID)
    rownames(estpers.upr) <- predvalues$time
    return(list(time=predvalues$time, persistence.prob=estpers, lower=estpers.lwr, upper=estpers.upr))
    }
  if(nlevels(dat$turbineID)==1){
    surv.mod<-coxph(Surv(perstime, status)~1, data=dat)
    predvalues <- summary(survfit(surv.mod))
    estpers <- predvalues$surv
    names(estpers) <- predvalues$time
    estpers.lwr <- predvalues$lower
    names(estpers.lwr) <- predvalues$time
    estpers.upr <- predvalues$upper
    names(estpers.upr) <- predvalues$time
    return(list(time=predvalues$time, persistence.prob=estpers, lower=estpers.lwr, upper=estpers.upr))
    }
 }
if(pers.const){
  if(nlevels(dat$turbineID)>1){
    mod<-survreg(Surv(perstime, status)~turbineID, data=dat, dist="exponential")
    estperstime <- predict(mod, newdata=newdat, type="response")
    # alternative way: S <- exp(-1/estperstime)
    b<-coef(mod)
    X<-model.matrix(~turbineID, data=newdat)
    bx<-X%*%b
    lambda<-exp(-bx)
    S<-exp(-lambda)
    sderrors<-summary(mod)$table[,2]
    Sran<-matrix(ncol=nlevels(dat$turbineID), nrow=R)
    for(i in 1:R){
    bran<-rnorm(nlevels(dat$turbineID), b, sderrors)
    branx<-X%*%bran
    lambda<-exp(-(branx))
    Sran[i,]<-exp(-lambda)
    }

    pers.prob.upper95 <- apply(Sran, MARGIN=2,quantile, prob=0.975)
    pers.prob.lower95 <- apply(Sran, MARGIN=2,quantile, prob=0.025)
    return(data.frame(turbineID=levels(dat$turbineID), persistence.prob=S, lower=pers.prob.lower95,
      upper=pers.prob.upper95, mean.persistence.time=estperstime))
      }
  if(nlevels(dat$turbineID)==1){
    mod<-survreg(Surv(perstime, status)~1, data=dat, dist="exponential")
    estperstime <- predict(mod, newdata=newdat, type="response")
    # alternative way: S <- exp(-1/estperstime)
    b<-coef(mod)
    X<-model.matrix(~1, data=newdat)
    bx<-X%*%b
    lambda<-exp(-bx)
    S<-exp(-lambda)
    sderrors<-summary(mod)$table[,2]
    Sran<-numeric(R)
    for(i in 1:R){
    bran<-rnorm(1, b, sderrors)
    branx<-X%*%bran
    lambda<-exp(-(branx))
    Sran[i]<-exp(-lambda)
    }

    pers.prob.upper95 <- quantile(Sran, prob=0.975)
    pers.prob.lower95 <- quantile(Sran, prob=0.025)
    return(data.frame(persistence.prob=S, lower=pers.prob.lower95,
      upper=pers.prob.upper95, mean.persistence.time=estperstime))
      }  
  
    }      
}

