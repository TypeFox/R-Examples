require("DoseFinding")
## effect curve estimate for linear type models!!

########################################################################
#### Testing function to generate doses and sample size allocs.
genDFdats <- function(model, argsMod, doses, n, sigma, mu = NULL){
  nD <- length(doses)
  dose <- sort(doses)
  if (length(n) == 1) n <- rep(n, nD)
  dose <- rep(dose,  n)
  args <- c(list(dose), argsMod)
  mu <- do.call(model, args)
  data.frame(dose = dose, resp = mu + rnorm(sum(n), sd = sigma))
}

getDosSampSiz <- function(){
  # generate dose levels
  mD <- runif(1, 0, 1500)
  nD <- max(rpois(1, 5), 4)
  p <- rgamma(nD, 3)
  p <- cumsum(p/sum(p))
  doses <- signif(c(0, mD*p), 3)

  # sample size allocations
  totSS <- rpois(1, rexp(1, 1/250))
  totSS <- max(totSS, 50)
  p <- rgamma(nD+1, 3);p <- p/sum(p)
  n <- round(p*totSS)
  n[n==0] <- rpois(sum(n==0), 1)+1
  list(doses=doses, n=n)
}

getDFdataSet <- function(doses, n){
  if(missing(doses) & missing(n)){
    ll <- getDosSampSiz()    
  } else {
    ll <- list(doses = doses, n=n)
  }

  e0 <- rnorm(1, 0, 10)
  eMax <- rgamma(1, abs(e0)*0.5, 0.5)
  sig <- eMax/runif(1, 0.5, 5)
  if(runif(1)<0.3){
    aa <- genDFdats("betaMod", c(e0 = e0, eMax = eMax, delta1=runif(1, 0.5, 4),
                delta2=runif(1, 0.5, 4), scal=1.2*max(ll$doses)),
                ll$doses, ll$n, sig)
  } else {
    aa <- genDFdats("sigEmax", c(e0 = e0, eMax = eMax,
                                 ed50=runif(1, 0.05*max(ll$doses), 1.5*max(ll$doses)),
                                 h=runif(1, 0.5, 4)), ll$doses, ll$n, sig)    
  }
  N <- sum(ll$n)
  center <- c("blue", "green", "red", "yellow", "silver")
  aa <- data.frame(x= aa$dose, y=aa$resp, center=as.factor(sample(center, N, replace = T)),
                   age=runif(N, 1, 100))
  aa[sample(1:nrow(aa)),]  
}

########################################################################
########################################################################
#### Generate data sets and compare results of fitDRModel
#### to the result of nls and lm for AIC function (if these are consistent
#### parameter estimates, residual sum of square and degrees of freedom are
#### consistent) and the vcov function (if these are consistent parameter
#### estimates, RSS, df and gradient are consistent)
########################################################################

########################################################################
#### beta Model
set.seed(2000)
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)

# without covariates
bnds <- matrix(c(0.05, 0.05, 6, 6), nrow=2)
fit0 <- fitMod(x, y, datset, model = "betaMod", addCovars = ~1,
               addArgs=list(scal=1.2*max(datset$x)), bnds=bnds, start=c(0.6, 0.6))
fitnls <- nls(y~betaMod(x, e0, emax, delta1, delta2, 1.2*max(datset$x)),
              start=c(e0=15, emax=14, delta1=0.8, delta2=0.5), data=datset)
AIC(fit0)
AIC(fitnls)
summary(fit0)
summary(fitnls)

vcov(fit0)
vcov(fitnls)

predict(fit0, predType="effect-curve", se.fit=TRUE)
predict(fit0, predType="full-model", se.fit=TRUE)

TD(fit0, Delta = 1)

# with covariates
fit0 <- fitMod(x, y, datset, model="betaMod", addCovars = ~age+center,
               addArgs=list(scal=1.2*max(datset$x)), bnds=bnds)
XX <- model.matrix(~center+age, data=datset)
scl <- 1.2*max(datset$x)
fitnls <- nls(y~cbind(XX, betaMod(x, 0, 1, delta1, delta2, scl)),
              data=datset, start=c(delta1=1, delta2=0.2),
              algorithm = "plinear")
AIC(fit0)
AIC(fitnls)
summary(fit0)
summary(fitnls)

vcov(fit0 )
vcov(fitnls)

predict(fit0, predType="effect-curve", doseSeq = c(0, 100), se.fit=T)

predict(fit0, predType="full-model", se.fit=T,
        newdata = data.frame(x = c(0,100), center = as.factor("yellow"), age = 50))

TD(fit0, Delta = 1)


########################################################################
#### emax Model
set.seed(15)
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)

# without covariates
bnds <- c(1e-5, max(datset$x))
fit0 <- fitMod(x,y, datset, model="emax", addCovars = ~1, bnds=bnds)
fitnls <- nls(y~emax(x, e0, emax, ed50),
              start=c(e0=-1, emax=1.3, ed50=0.1), data=datset)
AIC(fit0)
AIC(fitnls)
summary(fit0)
summary(fitnls)

vcov(fit0 )
vcov(fitnls)

predict(fit0, predType="effect-curve", se.fit=T)

predict(fit0, predType="full-model", se.fit=T)

TD(fit0, Delta = 0.005)

# with covariates
fit0 <- fitMod(x,y, datset, model="emax", addCovars = ~age+center, bnds=bnds)
XX <- model.matrix(~center+age, data=datset)
fitnls <- nls(y~cbind(XX, emax(x, 0, 1, ed50)),
              data=datset, start=list(ed50=1), algorithm = "plinear")
AIC(fit0)
AIC(fitnls)
summary(fit0)
summary(fitnls)

vcov(fit0 )
vcov(fitnls)

predict(fit0, predType="effect-curve", doseSeq = c(0, 100), se.fit=T)

predict(fit0, predType="full-model", se.fit=T,
        newdata = data.frame(x = c(0,100), center = as.factor("silver"), age = 50))

TD(fit0, Delta = 0.005)

########################################################################
#### sigEmax Model 
## set.seed(25) # example where nls and bndnls find different optimum
set.seed(13)
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)

# without covariates
bnds <- matrix(c(1e-5, 1e-5, max(datset$x), 30), nrow=2)
fit0 <- fitMod(x,y, datset, model = "sigEmax", addCovars = ~1, bnds=bnds)
fitnls <- nls(y~sigEmax(x, e0, emax, ed50, h),
              start=c(e0=6, emax=17, ed50=240, h=2), data=datset)
AIC(fit0)
AIC(fitnls)
summary(fit0)
summary(fitnls)

vcov(fit0 )
vcov(fitnls)

predict(fit0, predType="effect-curve", se.fit=T)

predict(fit0, predType="full-model", se.fit=T)

TD(fit0, Delta = 1)

# with covariates
fit0 <- fitMod(x,y, datset, model="sigEmax", addCovars = ~age+center, bnds=bnds)
XX <- model.matrix(~center+age, data=datset)
fitnls <- nls(y~cbind(XX, sigEmax(x, 0, 1, ed50, h)),
              data=datset, start=list(ed50=368, h=2), algorithm = "plinear")
AIC(fit0)
AIC(fitnls)
summary(fit0)
summary(fitnls)

vcov(fit0 )
vcov(fitnls)

predict(fit0, predType="effect-curve", doseSeq = c(0, 100), se.fit=T)

predict(fit0, predType="full-model", se.fit=T,
        newdata = data.frame(x = c(0,100), center = as.factor("silver"), age = 50))

TD(fit0, Delta = 1)

########################################################################
#### logistic Model
set.seed(200)
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)

# without covariates 
bnds <- matrix(c(1e-5, 1e-5, max(datset$x), max(datset$x)/2), nrow=2)
fit0 <- fitMod(x,y, datset, model="logistic", addCovars = ~1, bnds=bnds)
fitnls <- nls(y~logistic(x, e0, emax, ed50, delta),
              start=c(e0=0, emax=16, ed50=250, delta=90), data=datset)
AIC(fit0)
AIC(fitnls)
summary(fit0)
summary(fitnls)

vcov(fit0 )
vcov(fitnls)

predict(fit0, predType="effect-curve", se.fit=T)

predict(fit0, predType="full-model", se.fit=T)

TD(fit0, Delta = 0.5)

# with covariates (example where nls and bndnls find different optima)
fit0 <- fitMod(x,y, datset, model="logistic", addCovars = ~age+center, bnds=bnds) 
XX <- model.matrix(~center+age, data=datset)
fitnls <- nls(y~cbind(XX, logistic(x, 0, 1, ed50, delta)),
              data=datset, start=list(ed50=220, delta=48), algorithm = "plinear")
AIC(fit0)
AIC(fitnls)
summary(fit0)
summary(fitnls)

vcov(fit0 )
vcov(fitnls)

predict(fit0, predType="effect-curve", doseSeq = c(0, 100), se.fit=T)

predict(fit0, predType="full-model", se.fit=T,
        newdata = data.frame(x = c(0,100), center = as.factor("silver"), age = 5))

TD(fit0, Delta = 0.02)

########################################################################
#### exponential Model
set.seed(4)
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)

# without covariates
bnds <- c(0.1, 2)*max(datset$x)
fit0 <- fitMod(x,y, datset, model = "exponential", addCovars = ~1, bnds=bnds)
fitnls <- nls(y~exponential(x, e0, e1, delta),
              start=coef(fit0), data=datset)
AIC(fit0)
AIC(fitnls)
summary(fit0)
summary(fitnls)

vcov(fit0 )
vcov(fitnls)

predict(fit0, predType="effect-curve", se.fit=T)

predict(fit0, predType="full-model", se.fit=T)

TD(fit0, Delta = 0.1)

# with covariates
bnds <- c(0.1, 2)*max(datset$x)
fit0 <- fitMod(x,y, datset, model = "exponential", addCovars = ~age+center,
               bnds=bnds)
XX <- model.matrix(~center+age, data=datset)
fitnls <- nls(y~cbind(XX, exponential(x, 0, 1, delta)),
              data=datset, start=c(delta=450), algorithm = "plinear")
AIC(fit0)
AIC(fitnls)
summary(fit0)
summary(fitnls)

vcov(fit0 )
vcov(fitnls)

predict(fit0, predType="effect-curve", doseSeq = c(0, 100), se.fit=T)

predict(fit0, predType="full-model", se.fit=T,
        newdata = data.frame(x = c(0,100), center = as.factor("blue"), age = 50))

TD(fit0, Delta = 0.1)

########################################################################
#### linear model
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)

# without covariates
fit0 <- fitMod(x,y, datset, model = "linear", addCovars = ~1)
fitlm <- lm(y~x, data=datset)
AIC(fit0)
AIC(fitlm)
summary(fit0)
summary(fitlm)

vcov(fit0 )
vcov(fitlm)

predict(fit0, predType="effect-curve", se.fit=T)

TD(fit0, Delta = 1)

# with covariates
fit0 <- fitMod(x,y, datset, model = "linear", addCovars = ~age+center)
fitlm <- lm(y~x+age+center, data=datset)
AIC(fit0)
AIC(fitlm)
summary(fit0) 
summary(fitlm) 

vcov(fit0 )
vcov(fitlm)

predict(fit0, predType="effect-curve", se.fit=T)
predict(fit0, predType = "f", se.fit = T,
        newdata = data.frame(x=c(0,1,2,100), age = 30, center = as.factor("blue")))
predict(fitlm, se.fit = T,
        newdata = data.frame(x=c(0,1,2,100), age = 30, center = as.factor("blue")))

TD(fit0, Delta = 1)

########################################################################
#### linlog model
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)
off <- 0.05*max(datset$x)

# without covariates
fit0 <- fitMod(x,y, datset, model = "linlog", addCovars = ~1,addArgs=list(off=off))
fitlm <- lm(y~log(x+off), data=datset)
AIC(fit0)
AIC(fitlm)
summary(fit0)
summary(fitlm)

vcov(fit0 )
vcov(fitlm)

predict(fit0, predType="effect-curve", se.fit=T) ## bug ##

TD(fit0, Delta = 1)

# with covariates
fit0 <- fitMod(x,y, datset, model = "linlog", addCovars = ~age+center,
               addArgs=list(off=off))
fitlm <- lm(y~log(x+off)+age+center, data=datset)
AIC(fit0)
AIC(fitlm)
summary(fit0)
summary(fitlm)

vcov(fit0 )
vcov(fitlm)

predict(fit0, predType = "f", se.fit = T, ## degrees of freedom wrong ##
        newdata = data.frame(x=c(0,1,2,100), age = 35, center = as.factor("blue")))
predict(fitlm, se.fit = T,
        newdata = data.frame(x=c(0,1,2,100), age = 35, center = as.factor("blue")))

TD(fit0, Delta = 1)


########################################################################
#### quadratic model
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)

# without covariates
fit0 <- fitMod(x,y, datset, model = "quadratic", addCovars = ~1)
fitlm <- lm(y~x+I(x^2), data=datset)
AIC(fit0)
AIC(fitlm)
summary(fit0)
summary(fitlm)

vcov(fit0 )
vcov(fitlm)

predict(fit0, predType="effect-curve", se.fit=T)

predict(fit0, predType="full-model", se.fit=T,
        newdata=data.frame(x=c(0, 10, 100)))
predict(fitlm,  se.fit=T,
        newdata=data.frame(x=c(0, 10, 100)))

TD(fit0, Delta = 1)

# with covariates
fit0 <- fitMod(x,y, datset, model = "quadratic", addCovars = ~age+center)
fitlm <- lm(y~x+I(x^2)+age+center, data=datset)
AIC(fit0)
AIC(fitlm)
summary(fit0)
summary(fitlm)

vcov(fit0 )
vcov(fitlm)

predict(fit0, predType = "f", se.fit = T,
        newdata=data.frame(x=c(0, 10, 100), age = 30, center = as.factor("blue")))
predict(fitlm, se.fit = T,
        newdata=data.frame(x=c(0, 10, 100), age = 30, center = as.factor("blue")))

TD(fit0, Delta = 0.1)

########################################################################
## ensure that predict with no argument uses the original data not the
## sorted data that were used for fitting

data(IBScovars)
ff <- fitMod(dose, resp, data=IBScovars, model="quadratic",
             addCovars = ~gender)
## should be all zero
predict(ff, predType = "ls-means")-
predict(ff, predType = "ls-means", doseSeq = IBScovars[,3])
predict(ff, predType = "full-model")-
predict(ff, predType = "full-model", newdata = IBScovars[,-2])
predict(ff, predType = "effect-curve")-
predict(ff, predType = "effect-curve", doseSeq = IBScovars[,3])

ff2 <- fitMod(dose, resp, data=IBScovars, model="quadratic")
## should be all zero
predict(ff2, predType = "ls-means")-
predict(ff2, predType = "ls-means", doseSeq = IBScovars[,3])
predict(ff2, predType = "full-model")-
predict(ff2, predType = "full-model", newdata = IBScovars[,-2])
predict(ff2, predType = "effect-curve")-
predict(ff2, predType = "effect-curve", doseSeq = IBScovars[,3])

dose <- unique(IBScovars$dose)
ord <- c(2,4,1,3,5)
mns <- tapply(IBScovars$resp, IBScovars$dose, mean)[ord]
ff3 <- fitMod(dose, mns, S=diag(5), model="quadratic", type = "general")
predict(ff3, predType = "ls-means")-
predict(ff3, predType = "ls-means", doseSeq = dose)
predict(ff3, predType = "effect-curve")-
predict(ff3, predType = "effect-curve", doseSeq = dose)

########################################################################
## ensure that S is also sorted when the dose is not entered sorted
dose <- sort(unique(IBScovars$dose))
mns <- tapply(IBScovars$resp, IBScovars$dose, mean)
S <- c(1000,1,1,1,1)*diag(5)
ff1 <- fitMod(dose, mns, S = S, model="linear", type="general")
## fit unsorted
dose <- unique(IBScovars$dose)
ord <- c(2,4,1,3,5)
mns <- tapply(IBScovars$resp, IBScovars$dose, mean)[ord]
ff2 <- fitMod(dose, mns, S = S, model="linear", type="general")
ff3 <- fitMod(dose, mns, S = S[ord,ord], model="linear", type="general")
## coef(ff1) & coef(ff3) should be equal
coef(ff1)
coef(ff3)
