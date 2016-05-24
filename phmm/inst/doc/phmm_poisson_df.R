### R code from vignette source 'phmm_poisson_df.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: phmm_poisson_df.Rnw:36-62
###################################################
options(width=75)
library(phmm)

n <- 50      # total sample size
nclust <- 5  # number of clusters
clusters <- rep(1:nclust,each=n/nclust)
beta0 <- c(1,2)
set.seed(13)

Z <-cbind(Z1=sample(0:1,n,replace=TRUE),
          Z2=sample(0:1,n,replace=TRUE),
          Z3=sample(0:1,n,replace=TRUE))
b <- cbind(rep(rnorm(nclust), each=n/nclust),
           rep(rnorm(nclust), each=n/nclust))
Wb <- matrix(0,n,2)
for( j in 1:2) Wb[,j] <- Z[,j]*b[,j]
Wb <- apply(Wb,1,sum)
T <- -log(runif(n,0,1))*exp(-Z[,c('Z1','Z2')]%*%beta0-Wb)
C <- runif(n,0,1)
time <- ifelse(T<C,T,C)
event <- ifelse(T <= C,1,0)
sum(event)
phmmd <- data.frame(Z)
phmmd$cluster <- clusters
phmmd$time <- time
phmmd$event <- event


###################################################
### code chunk number 2: phmm_poisson_df.Rnw:66-71
###################################################
fit.ph <- coxph(Surv(time, event) ~ Z1 + Z2, 
   phmmd, method="breslow", x=TRUE, y=TRUE)

summary(fit.ph)
fit.ph$loglik[2]


###################################################
### code chunk number 3: phmm_poisson_df.Rnw:88-101
###################################################
ppd <- as.data.frame(as.matrix(pseudoPoisPHMM(fit.ph)))

# pois likelihood
poisl <- c()
eventtimes <- sort(phmmd$time[phmmd$event == 1])
for(h in 1:length(eventtimes)){
  js <- ppd$time == eventtimes[h] & ppd$m >= 1  # j star
  j  <- ppd$time == eventtimes[h]
  if(sum(js) > 1) stop("tied event times")
  poisl <- c(poisl, 
    ppd[js, "N"]*exp(-1)*exp(ppd[js, "linear.predictors"])/
    sum(ppd[j, "N"]*exp(ppd[j, "linear.predictors"])))
}


###################################################
### code chunk number 4: phmm_poisson_df.Rnw:104-107
###################################################
sum(log(poisl))

sum(log(poisl)) - fit.ph$loglik[2]


###################################################
### code chunk number 5: phmm_poisson_df.Rnw:110-111
###################################################
length(fit.ph$coef) + sum(phmmd$event)


###################################################
### code chunk number 6: phmm_poisson_df.Rnw:116-124
###################################################
ppd$t <- as.factor(ppd$time)
fit.glm <- glm(m~-1+t+z1+z2+offset(log(N)), 
  ppd, family=poisson)

summary(fit.glm)
fit.ph$coef
logLik(fit.glm)
logLik(fit.glm)[1] - sum(log(poisl))


###################################################
### code chunk number 7: phmm_poisson_df.Rnw:128-130
###################################################
bh <- basehaz(fit.ph, centered = FALSE)
log(bh$hazard - c(0,bh$hazard[1:(length(bh$hazard)-1)]))[1:10]


###################################################
### code chunk number 8: phmm_poisson_df.Rnw:135-139
###################################################
fit.phmm <- phmm(Surv(time, event) ~ Z1 + Z2 + (Z1 + Z2|cluster), 
   phmmd, Gbs = 100, Gbsvar = 1000, VARSTART = 1,
   NINIT = 10, MAXSTEP = 100, CONVERG=90)
summary(fit.phmm)


###################################################
### code chunk number 9: phmm_poisson_df.Rnw:143-155
###################################################
ppd <- as.data.frame(as.matrix(pseudoPoisPHMM(fit.phmm)))

poisl <- c()
eventtimes <- sort(phmmd$time[phmmd$event == 1])
for(h in 1:length(eventtimes)){
  js <- ppd$time == eventtimes[h] & ppd$m >= 1  # j star
  j  <- ppd$time == eventtimes[h]
  if(sum(js) > 1) stop("tied event times")
  poisl <- c(poisl, 
    ppd[js, "N"]*exp(-1)*exp(ppd[js, "linear.predictors"])/
    sum(ppd[j, "N"]*exp(ppd[j, "linear.predictors"])))
}


###################################################
### code chunk number 10: phmm_poisson_df.Rnw:158-161
###################################################
sum(log(poisl))

sum(log(poisl)) - fit.phmm$loglik[1]


###################################################
### code chunk number 11: phmm_poisson_df.Rnw:164-166
###################################################
# Poisson GLMM degrees of freedom  length(unique(x$cluster)) * x$nrandom + x$nfixed
traceHat(fit.phmm, "pseudoPois") # + 2*sum(phmmd$event)


###################################################
### code chunk number 12: phmm_poisson_df.Rnw:171-182
###################################################
library(lme4)
ppd$t <- as.factor(ppd$time)
fit.lmer <- lmer(m~-1+t+z1+z2+
  (z1+z2|cluster)+offset(log(N)), 
  data=ppd, family=poisson)

summary(fit.lmer)$coef
fit.phmm$coef
logLik(fit.lmer)

sum(log(poisl)) - logLik(fit.lmer)[1]


###################################################
### code chunk number 13: phmm_poisson_df.Rnw:185-186
###################################################
log(fit.phmm$lambda)[1:10]


