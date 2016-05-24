### family-specific functions

# Here, nlog.double.eps=-log(.Machine[["double.eps"]]) defined in .onLoad is used to
# decide when exp(eta)==exp(eta)+1

## Bernoulli logit
lpY.Bernoulli.logit<-function(Y,eta,dispersion=NULL,fam.par=NULL){
  ifelse(eta>=.latentnetEnv$nlog.double.eps,eta*(Y-1),eta*Y-log1p(exp(eta)))
}
lpYc.Bernoulli.logit<-function(Y,eta,dispersion=NULL,fam.par=NULL){
  ifelse(eta>=.latentnetEnv$nlog.double.eps,eta*(Y-1),eta*Y-log1p(exp(eta)))
}
pY.Bernoulli.logit<-function(Y=1,eta,dispersion=NULL,fam.par=NULL){
  ifelse(eta>=.latentnetEnv$nlog.double.eps,exp(eta*(Y-1)),exp(eta*Y)/(exp(eta)+1))
}
dlpY.deta.Bernoulli.logit<-function(Y,eta,dispersion=NULL,fam.par=NULL) Y-EY.Bernoulli.logit(eta,fam.par)
rsm.Bernoulli.logit<-function(eta,dispersion=NULL,fam.par=NULL){
  n<-dim(eta)[1]
  matrix(rbinom(n*n,1,pY.Bernoulli.logit(eta=eta,dispersion=NULL,fam.par=NULL)),n,n)
}
EY.Bernoulli.logit<-function(eta,dispersion=NULL,fam.par=NULL) 1/(1+exp(-eta))

## Binomial logit

lpY.binomial.logit<-function(Y,eta,dispersion=NULL,fam.par){
  ifelse(eta>=.latentnetEnv$nlog.double.eps,eta*(Y-fam.par[["trials"]])+lchoose(fam.par[["trials"]],Y),
         eta*Y-fam.par[["trials"]]*log1p(exp(eta))+lchoose(fam.par[["trials"]],Y))
}
lpYc.binomial.logit<-function(Y,eta,dispersion=NULL,fam.par){
  ifelse(eta>=.latentnetEnv$nlog.double.eps,eta*(Y-fam.par[["trials"]]),
         eta*Y-fam.par[["trials"]]*log1p(exp(eta)))
}
pY.binomial.logit<-function(Y,eta,dispersion=NULL,fam.par) exp(lpY.binomial.logit(Y,eta,dispersion=NULL,fam.par))
dlpY.deta.binomial.logit<-function(Y,eta,dispersion=NULL,fam.par) (Y-EY.binomial.logit(eta,dispersion=NULL,fam.par))
rsm.binomial.logit<-function(eta,dispersion=NULL,fam.par){
  n<-dim(eta)[1]
  matrix(rbinom(n*n,fam.par[["trials"]],EY.binomial.logit(eta,dispersion=NULL,fam.par)/fam.par[["trials"]]),n,n)
}
EY.binomial.logit<-function(eta,dispersion=NULL,fam.par) fam.par[["trials"]]/(1+exp(-eta))

## Poisson log

lpY.Poisson.log<-function(Y,eta,dispersion=NULL,fam.par=NULL) dpois(Y,EY.Poisson.log(eta,dispersion=NULL,fam.par),TRUE)
lpYc.Poisson.log<-function(Y,eta,dispersion=NULL,fam.par=NULL) Y*eta-exp(eta)
pY.Poisson.log<-function(Y,eta,dispersion=NULL,fam.par=NULL) dpois(Y,EY.Poisson.log(eta,dispersion=NULL,fam.par),FALSE)
dlpY.deta.Poisson.log<-function(Y,eta,dispersion=NULL,fam.par=NULL) Y-EY.Poisson.log(eta,dispersion=NULL,fam.par)
rsm.Poisson.log<-function(eta,dispersion=NULL,fam.par=NULL){
  n<-dim(eta)[1]
  matrix(rpois(n*n,exp(eta)),n,n)
}
EY.Poisson.log<-function(eta,dispersion=NULL,fam.par=NULL) exp(eta)

## normal identity

lpY.normal.identity<-function(Y,eta,dispersion=NULL,fam.par=NULL) dnorm(Y,eta,sqrt(dispersion),TRUE)
lpYc.normal.identity<-function(Y,eta,dispersion=NULL,fam.par=NULL) -(Y-eta)^2/dispersion/2-log(dispersion)/2
pY.normal.identity<-function(Y,eta,dispersion=NULL,fam.par=NULL) dnorm(Y,eta,sqrt(dispersion),FALSE)
dlpY.deta.normal.identity<-function(Y,eta,dispersion=NULL,fam.par=NULL) (Y-eta)/dispersion
dlpY.ddispersion.normal.identity<-function(Y,eta,dispersion=NULL,fam.par=NULL) sum(na.omit(c(Y-eta))^2/dispersion^2/2-1/dispersion/2)
rsm.normal.identity<-function(eta,dispersion=NULL,fam.par=NULL){
  n<-dim(eta)[1]
  matrix(rnorm(n*n,eta,sqrt(dispersion)),n,n)
}
EY.normal.identity<-function(eta,dispersion=NULL,fam.par=NULL) eta

## Dispatcher functions

family.IDs<-list(Bernoulli=1,
                 Bernoulli.logit=1,
                 binomial=2,
                 binomial.logit=2,
                 Poisson=3,
                 Poisson.log=3,
                 Bernoulli.cont.logit=4,
                 binomial.cont.logit=5,
                 Poisson.cont.log=6,
                 normal=7,
                 normal.identity=7,
                 Gaussian=7,
                 Gaussian.identity=7)

family.names<-c("Bernoulli.logit",
                "binomial.logit",
                "Poisson.log",
                "Bernoulli.cont.logit",
                "binomial.cont.logit",
                "Poisson.cont.log",
                "normal.identity")


lpY.fs<-c(lpY.Bernoulli.logit,
          lpY.binomial.logit,
          lpY.Poisson.log,
          lpY.Bernoulli.logit,
          lpY.binomial.logit,
          lpY.Poisson.log,
          lpY.normal.identity)

lpYc.fs<-c(lpYc.Bernoulli.logit,
           lpYc.binomial.logit,
           lpYc.Poisson.log,
           lpYc.Bernoulli.logit,
           lpYc.binomial.logit,
           lpYc.Poisson.log,
           lpYc.normal.identity)
          
pY.fs<-c(pY.Bernoulli.logit,
         pY.binomial.logit,
         pY.Poisson.log,
         pY.Bernoulli.logit,
         pY.binomial.logit,
         pY.Poisson.log,
         pY.normal.identity)

dlpY.deta.fs<-c(dlpY.deta.Bernoulli.logit,
                dlpY.deta.binomial.logit,
                dlpY.deta.Poisson.log,
                dlpY.deta.Bernoulli.logit,
                dlpY.deta.binomial.logit,
                dlpY.deta.Poisson.log,
                dlpY.deta.normal.identity)

dlpY.ddispersion.fs<-c(function(...) 0,
                       function(...) 0,
                       function(...) 0,
                       function(...) 0,
                       function(...) 0,
                       function(...) 0,
                       dlpY.ddispersion.normal.identity)

rsm.fs<-c(rsm.Bernoulli.logit,
          rsm.binomial.logit,
          rsm.Poisson.log,
          rsm.Bernoulli.logit,
          rsm.binomial.logit,
          rsm.Poisson.log,
          rsm.normal.identity)
          
EY.fs<-c(EY.Bernoulli.logit,
         EY.binomial.logit,
         EY.Poisson.log,
         EY.Bernoulli.logit,
         EY.binomial.logit,
         EY.Poisson.log,
         EY.normal.identity)

fam.par.check<-function(model){
  if(model[["familyID"]]==2){
    if(is.null(model[["fam.par"]][["trials"]]))
      stop("Binomial family requires parameter `trials'.")
    model[["iconsts"]]<-model[["fam.par"]][["trials"]]
  }
  if(model[["familyID"]]==5){
    if(is.null(model[["fam.par"]][["trials"]]))
      stop("Binomial (cont) family requires parameter `trials'.")
    model[["dconsts"]]<-model[["fam.par"]][["trials"]]
  }
  if(model[["familyID"]]==7){
    if(is.null(model[["fam.par"]][["prior.var"]])||is.null(model[["fam.par"]][["prior.var.df"]]))
      stop("Normal family requires prior parameters `prior.var` and `prior.var.df`.")
    model[["prior"]][["dispersion"]]<-model[["fam.par"]][["prior.var"]]
    model[["prior"]][["dispersion.df"]]<-model[["fam.par"]][["prior.var.df"]]
    model[["dispersion"]]<-TRUE
  }
  model
}
