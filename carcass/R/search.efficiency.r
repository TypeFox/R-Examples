###########################################################################
# function to estimate searcher efficiency with 95% confidence interval 
# based on experimental data
###########################################################################

# Protocol of changes
# - fk, 2.7.14:  line 152: change npersons==2 into npersons<=4 & npersons>1 (bug reported by Joanna Bernardino)
# - fk, 25.10.13: line 55: b changed to bsim@fixef[j,], so that uncertainty of estimates for fixed effects is also considered 
# - tr, 13.06.13: Data can either be provided with seperate vectors or as a data.frame

search.efficiency <- function(dat=NA, person=NA, visibility=NA, detected=NA, notdetected=NA, nsim=1000){
# A data.frame containing the following variables: 
  #$ person     : names of the persons who searched
  #$ visibility : visibility class
  #$ detected   : number of detected items
  #$ notdetected: number of not detected items
# Or alternatively, the following vectors:
# - person     : names of the persons who searched
# - visibility : visibility class
# - detected   : number of detected items
# - notdetected: number of not detected items
# Finally:
# - nsim: number of simulations to be drawn from the posterior distributions to describe the 95% credible intervals
  
if(is.na(person[1]) & class(dat) != "data.frame") {
  stop("Please provide the data either as a data.frame containing all data or, alternatively
as seperate vectors (i.e. the vector person, visibility, detected, nondetected).")
}
    
if(class(dat) != "data.frame") {
  if(length(person)!=length(visibility) | length(person)!=length(detected) | length(person)!=length(notdetected)) {
    stop("The vector 'person', 'visibility', 'detected' and 'notdetected' should all have the same length.")
  }
  dat <- data.frame(list(person=person, visibility=visibility, detected=detected, notdetected=notdetected))
}

dat$visibility <- factor(dat$visibility, levels=levels(dat$visibility))[drop=TRUE]   
dat$person <- factor(dat$person)[drop=TRUE]

npersons <- nlevels(dat$person)
nvisclass <- nlevels(dat$visibility)

if(npersons>2 & nvisclass>1){
  mod <- glmer(cbind(detected, notdetected) ~ visibility + (1|person), data=dat, family=binomial)
  newdat <- expand.grid(visibility=factor(levels(dat$visibility), levels=levels(dat$visibility)), person=levels(dat$person))
  b <- fixef(mod)
  bperson <- matrix(b, nrow=npersons, ncol=nvisclass, byrow=TRUE)
  bperson[,1] <- bperson[,1] + ranef(mod)$person[,1]
  rownames(bperson) <- levels(dat$person)
  newdat$f <- NA
  Xmat <- model.matrix(~visibility, data=newdat[newdat$person==levels(newdat$person)[1],]) 
  for(i in levels(dat$person)) newdat$f[newdat$person==i] <- plogis(Xmat%*%as.numeric(bperson[i,]))
   
  bsim <- sim(mod, n.sim=nsim)
  
  predmat <- matrix(nrow=nrow(newdat), ncol=nsim)
  for(i in levels(dat$person)){ 
    for(j in 1:nsim) {
    bpersi <- bsim@fixef[j,] 
    bpersi[1] <- bpersi[1] + bsim@ranef$person[j,i,1]
    predmat[newdat$person==i, j] <- plogis(Xmat%*%bpersi)
    }
    }
  newdat$lwr <- apply(predmat, 1, quantile, prob=0.025)
  newdat$upr <- apply(predmat, 1, quantile, prob=0.975)
  newdat$se <- apply(predmat, 1, sd)
  
  overall <- data.frame(visibility=factor(levels(dat$visibility), levels=levels(dat$visibility)))
  Xmat <- model.matrix(~visibility, data=overall)
  overall$f <- plogis(Xmat%*%b)
  predmat <- matrix(nrow=nrow(overall), ncol=nsim)
  for(i in 1:nsim) predmat[,i] <- plogis(Xmat%*%bsim@fixef[i,])
  overall$lwr <- apply(predmat, 1, quantile, prob=0.025)
  overall$upr <- apply(predmat, 1, quantile, prob=0.975)
  overall$se <- apply(predmat, 1, sd)
 
  }



if(npersons>1 & npersons<5 & nvisclass>1){
  mod <- glm(cbind(detected, notdetected) ~ visibility + person, data=dat, family=binomial)
  newdat <- expand.grid(visibility=factor(levels(dat$visibility), levels=levels(dat$visibility)), person=levels(dat$person))
  b <- coef(mod)
  Xmat <- model.matrix(~visibility + person, data=newdat) 
  newdat$f <- plogis(Xmat%*%b)
   
  bsim <- sim(mod, n.sim=nsim)
  
  predmat <- matrix(nrow=nrow(newdat), ncol=nsim)
  for(i in 1:nsim) predmat[, i] <- plogis(Xmat%*%bsim@coef[i,])

  newdat$lwr <- apply(predmat, 1, quantile, prob=0.025)
  newdat$upr <- apply(predmat, 1, quantile, prob=0.975)
  newdat$se <- apply(predmat, 1, sd)
  
  overall <- data.frame(visibility=factor(levels(dat$visibility), levels=levels(dat$visibility)))
  predmat1 <- matrix(nrow=nrow(overall), ncol=nsim)
  for(i in 1:nsim) predmat1[,i] <- tapply(predmat[,i], newdat$visibility, mean)
  overall$f <- apply(predmat1, 1, mean)
  overall$lwr <- apply(predmat1, 1, quantile, prob=0.025)
  overall$upr <- apply(predmat1, 1, quantile, prob=0.975)
  overall$se <- apply(predmat1, 1, sd)
 
  }
     

if(npersons==1 & nvisclass>1){
  mod <- glm(cbind(detected, notdetected) ~ visibility, data=dat, family=binomial)
  newdat <- expand.grid(visibility=factor(levels(dat$visibility), levels=levels(dat$visibility)), person=levels(dat$person))
  b <- coef(mod)
  Xmat <- model.matrix(~visibility, data=newdat) 
  newdat$f <- plogis(Xmat%*%b)
   
  bsim <- sim(mod, n.sim=nsim)
  
  predmat <- matrix(nrow=nrow(newdat), ncol=nsim)
  for(i in 1:nsim) predmat[, i] <- plogis(Xmat%*%bsim@coef[i,])

  newdat$lwr <- apply(predmat, 1, quantile, prob=0.025)
  newdat$upr <- apply(predmat, 1, quantile, prob=0.975)
  newdat$se <- apply(predmat, 1, sd)
  
  overall <- newdat
  }




if(npersons>=5 & nvisclass==1){
  mod <- glmer(cbind(detected, notdetected) ~ 1 + (1|person), data=dat, family=binomial)
  newdat <- expand.grid(visibility=levels(dat$visibility), person=levels(dat$person))
  b <- fixef(mod)
  bperson <- b + ranef(mod)$person
  newdat$f <- NA
  Xmat <- matrix(1, ncol=1, nrow=1) 
  for(i in levels(dat$person)) newdat$f[newdat$person==i] <- plogis(Xmat%*%as.numeric(bperson[i,]))
   
  bsim <- sim(mod, n.sim=nsim)
  
  predmat <- matrix(nrow=nrow(newdat), ncol=nsim)
  for(i in levels(dat$person)){ 
    for(j in 1:nsim) {
    bpersi <- bsim@fixef[j,] + bsim@ranef$person[j,i,1]
    predmat[newdat$person==i, j] <- plogis(Xmat%*%bpersi)
    }
    }
  newdat$lwr <- apply(predmat, 1, quantile, prob=0.025)
  newdat$upr <- apply(predmat, 1, quantile, prob=0.975)
  newdat$se <- apply(predmat, 1, sd)
  
  overall <- data.frame(visibility=levels(dat$visibility))
  Xmat <- model.matrix(~1, data=overall)
  overall$f <- plogis(Xmat%*%b)
  predmat <- matrix(nrow=nrow(overall), ncol=nsim)
  for(i in 1:nsim) predmat[,i] <- plogis(Xmat%*%bsim@fixef[i,])
  overall$lwr <- apply(predmat, 1, quantile, prob=0.025)
  overall$upr <- apply(predmat, 1, quantile, prob=0.975)
  overall$se <- apply(predmat, 1, sd)
  }



if(npersons<=4 & npersons>1 &nvisclass==1){
  mod <- glm(cbind(detected, notdetected) ~ person, data=dat, family=binomial)
  newdat <- expand.grid(visibility=levels(dat$visibility), person=levels(dat$person))
  b <- coef(mod)
  Xmat <- model.matrix(~1 + person, data=newdat) 
  newdat$f <- plogis(Xmat%*%b)
   
  bsim <- sim(mod, n.sim=nsim)
  
  predmat <- matrix(nrow=nrow(newdat), ncol=nsim)
  for(i in 1:nsim) predmat[, i] <- plogis(Xmat%*%bsim@coef[i,])

  newdat$lwr <- apply(predmat, 1, quantile, prob=0.025)
  newdat$upr <- apply(predmat, 1, quantile, prob=0.975)
  newdat$se <- apply(predmat, 1, sd)
  
  overall <- data.frame(visibility=levels(dat$visibility))
  predmat1 <- matrix(nrow=nrow(overall), ncol=nsim)
  for(i in 1:nsim) predmat1[,i] <- tapply(predmat[,i], newdat$visibility, mean)
  overall$f <- apply(predmat1, 1, mean)
  overall$lwr <- apply(predmat1, 1, quantile, prob=0.025)
  overall$upr <- apply(predmat1, 1, quantile, prob=0.975)
  overall$se <- apply(predmat1, 1, sd)
  }
     

if(npersons==1 & nvisclass==1){
  mod <- glm(cbind(detected, notdetected) ~ 1, data=dat, family=binomial)
  newdat <- expand.grid(visibility=levels(dat$visibility), person=levels(dat$person))
  b <- coef(mod)
  Xmat <- model.matrix(~1, data=newdat) 
  newdat$f <- plogis(Xmat%*%b)
   
  bsim <- sim(mod, n.sim=nsim)
  
  predmat <- matrix(nrow=nrow(newdat), ncol=nsim)
  for(i in 1:nsim) predmat[, i] <- plogis(Xmat%*%bsim@coef[i,])

  newdat$lwr <- apply(predmat, 1, quantile, prob=0.025)
  newdat$upr <- apply(predmat, 1, quantile, prob=0.975)
  newdat$se <- apply(predmat, 1, sd)
  
  overall <- newdat
  }

  return(list(f.perperson=newdat, f.average=overall))
  }
  


