require("DoseFinding")
if(!require("multcomp"))
  stop("need multcomp package to run this test")

########################################################################
#### multContTest
# functions to sample random DF data
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
  ll <- getDosSampSiz()
  e0 <- rnorm(1, 0, 10)
  eMax <- rgamma(1, abs(e0)*0.5, 0.5)*I(runif(1)<0.25)
  if(eMax > 0){ sig <- eMax/runif(1, 0.5, 5)}
  else { sig <- rgamma(1, abs(e0)*0.5, 0.5) }
  dosVec <- rep(ll$doses, ll$n)
  if(runif(1)<0.3){
    mnVec <- betaMod(dosVec, e0=e0, eMax=eMax, delta1=runif(1, 0.5, 5),
                     delta2=runif(1, 0.5, 5), scal=1.2*max(ll$doses))
  } else {
    mnVec <- logistic(dosVec, e0 = e0, eMax = eMax,
                      ed50=runif(1, 0.05*max(ll$doses), 1.5*max(ll$doses)),
                      delta=runif(1, 0.5, max(ll$doses)/2))
  }
  resp <- rnorm(sum(ll$n), mnVec, sig)
  N <- sum(ll$n)
  cov1 <- as.factor(rpois(N, 5))
  cov2 <- runif(N, 1, 100)
  aa <- data.frame(x= dosVec, y=resp, cov1=cov1, cov2=cov2)
  aa[sample(1:nrow(aa)),]  
}

#### simulate data and compare to output of glht of multcomp package and oldMCPMod function
set.seed(10)
dd <- getDFdataSet()
bet <- guesst(0.9*max(dd$x), p=0.8, "betaMod", scal = 1.2*max(dd$x),
              dMax = 0.7*max(dd$x), Maxd = max(dd$x))
sE <- guesst(c(0.5*max(dd$x), 0.7*max(dd$x)) , p=c(0.5, 0.9), "sigEmax")
models <- Mods(linear = NULL, betaMod = bet, sigEmax = sE,
               doses = sort(unique(dd$x)),
               addArgs=list(scal = 1.2*max(dd$x)))
obj <- MCTtest(x,y, dd, models=models, addCovars = ~cov1+cov2, pVal = T)
dd2 <- dd;dd2$x <- as.factor(dd$x)
fit <- lm(y~x+cov1+cov2, data=dd2)
mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
summary(mcp)
print(obj, digits = 3)

obj <- MCTtest(x,y, dd, models=models, addCovars = ~1, pVal = T)
dd2 <- dd;dd2$x <- as.factor(dd$x)
fit <- lm(y~x, data=dd2)
mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
summary(mcp)
print(obj, digits = 3)

#### different model set
set.seed(10)
dd <- getDFdataSet()
mD <- max(dd$x)
lg1 <- guesst(c(0.3*mD, 0.4*mD), c(0.3, 0.9), "logistic")
lg2 <- guesst(c(0.3*mD, 0.4*mD), c(0.3, 0.5), "logistic")
expo <- guesst(c(0.9*mD), c(0.7), "exponential", Maxd=mD)
quad <- guesst(c(0.6*mD), c(1), "quadratic")
models <- Mods(linlog = NULL, logistic = rbind(lg1, lg2),
                   exponential = expo, quadratic = quad,
                   doses = sort(unique(dd$x)), addArgs=list(off = 0.2*max(dd$x)))

obj <- MCTtest(x,y, dd, models=models, addCovars = ~cov1+cov2, pVal = T)
dd2 <- dd;dd2$x <- as.factor(dd$x)
fit <- lm(y~x+cov1+cov2, data=dd2)
mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
summary(mcp)
print(obj, digits = 3)

obj <- MCTtest(x,y, dd, models=models, addCovars = ~1, pVal = T)
dd2 <- dd;dd2$x <- as.factor(dd$x)
fit <- lm(y~x, data=dd2)
mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
summary(mcp)
print(obj, digits = 3)

#### contrast matrix handed over
set.seed(23)
dd <- getDFdataSet()
mD <- max(dd$x)
lg1 <- guesst(c(0.3*mD, 0.4*mD), c(0.3, 0.9), "logistic")
lg2 <- guesst(c(0.3*mD, 0.4*mD), c(0.3, 0.5), "logistic")
expo <- guesst(c(0.9*mD), c(0.7), "exponential", Maxd=mD)
quad <- guesst(c(0.6*mD), c(1), "quadratic")
models <- Mods(linlog = NULL, logistic = rbind(lg1, lg2),
                   exponential = expo, quadratic = quad,
                   doses = dd$x, addArgs=list(off = 0.2*max(dd$x)))

obj <- MCTtest(x,y, dd, models=models, addCovars = ~cov1+cov2, pVal = T)
contMat <- obj$contMat
obj <- MCTtest(x,y, dd, models=models, addCovars = ~1, pVal = T, contMat = contMat)
dd2 <- dd
dd2$x <- as.factor(dd2$x)
fit <- lm(y~x, data=dd2)
mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
summary(mcp)
obj

########################################################################
#### some binary test cases
getDFdataSet.bin <- function(doses, n){
  ll <- getDosSampSiz()
  ll$n <- ll$n+10
  e0 <- rnorm(1, 0, sqrt(3.28))
  eMax <- rnorm(1, 0, 5)
  dosVec <- rep(ll$doses, ll$n)
  if(runif(1)<0.3){
    mn <- betaMod(dosVec, e0 = e0, eMax = eMax, delta1=runif(1, 0.5, 5),
                  delta2=runif(1, 0.5, 5), scal=1.2*max(ll$doses))
  } else {
    mn <- logistic(dosVec, e0 = e0,
                   eMax = eMax, ed50=runif(1, 0.05*max(ll$doses), 1.5*max(ll$doses)),
                   delta=runif(1, 0.5, max(ll$doses)/2))
  }
  resp <- rbinom(length(ll$n), ll$n, 1/(1+exp(-mn)))
  aa <- data.frame(dose = ll$doses, resp = resp)    
  aa <- data.frame(x= aa$dose, y=aa$resp/ll$n, n=ll$n)
  aa[sample(1:nrow(aa)),]  
}

set.seed(1909)
dd <- getDFdataSet.bin()
bet <- guesst(0.9*max(dd$x), p=0.8, "betaMod", scal = 1.2*max(dd$x), dMax = 0.7*max(dd$x),
              Maxd = max(dd$x))
sE <- guesst(c(0.5*max(dd$x), 0.7*max(dd$x)) , p=c(0.5, 0.9), "sigEmax")
models <- Mods(linear = NULL, betaMod = bet, sigEmax = sE,
                   doses = sort(unique(dd$x)), addArgs=list(scal = 1.2*max(dd$x)))
logReg <- glm(y~as.factor(x)-1, family=binomial, data=dd, weights = n)
dePar <- coef(logReg)
vCov <- vcov(logReg)
dose <- sort(unique(dd$x))
obj <- MCTtest(dose, dePar, S=vCov, models=models, type="general",
               df=Inf, pVal = T)
dd2 <- dd;dd2$x <- as.factor(dd$x)
fit <- glm(y~x-1, family = binomial, data=dd2, weights = n)
mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
summary(mcp)
print(obj, digits = 3)

set.seed(1997)
dd <- getDFdataSet.bin()
bet <- guesst(0.9*max(dd$x), p=0.8, "betaMod", scal = 1.2*max(dd$x),
              dMax = 0.7*max(dd$x), Maxd = max(dd$x))
sE <- guesst(c(0.5*max(dd$x), 0.7*max(dd$x)) , p=c(0.5, 0.9), "sigEmax")
models <- Mods(linear = NULL, betaMod = bet, sigEmax = sE,direction = "decreasing",
                   addArgs=list(scal = 1.2*max(dd$x)), doses = sort(unique(dd$x)))
logReg <- glm(y~as.factor(x)-1, family=binomial, data=dd, weights = n)
dePar <- coef(logReg)
vCov <- vcov(logReg)
dose <- sort(unique(dd$x))
obj <- MCTtest(dose, dePar, S=vCov, models=models, type = "general",
               pVal = TRUE, df=Inf)
dd2 <- dd;dd2$x <- as.factor(dd$x)
fit <- glm(y~x-1, family = binomial, data=dd2, weights = n)
mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
summary(mcp)
print(obj, digits = 3)

set.seed(1)
dd <- getDFdataSet.bin()
bet <- guesst(0.9*max(dd$x), p=0.8, "betaMod", scal = 1.2*max(dd$x),
              dMax = 0.7*max(dd$x), Maxd = max(dd$x))
sE <- guesst(c(0.5*max(dd$x), 0.7*max(dd$x)) , p=c(0.5, 0.9), "sigEmax")
models <- Mods(linear = NULL, betaMod = bet, sigEmax = sE,
                   doses = sort(unique(dd$x)), addArgs=list(scal = 1.2*max(dd$x)))
logReg <- glm(y~as.factor(x)-1, family=binomial, data=dd, weights = n)
dePar <- coef(logReg)
vCov <- vcov(logReg)
dose <- sort(unique(dd$x))
obj <- MCTtest(dose, dePar, S=vCov, models=models, type = "general",
               pVal = T, df=Inf)
dd2 <- dd;dd2$x <- as.factor(dd$x)
fit <- glm(y~x-1, family = binomial, data=dd2, weights = n)
mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
summary(mcp)
print(obj, digits = 3)

## one-dimensional test
set.seed(1)
dd <- getDFdataSet.bin()
model <- Mods(linear = NULL, doses=sort(unique(dd$x)), addArgs=list(scal = 1.2*max(dd$x)))
logReg <- glm(y~as.factor(x)-1, family=binomial, data=dd, weights = n)
dePar <- coef(logReg)
vCov <- vcov(logReg)
dose <- sort(unique(dd$x))
obj <- MCTtest(dose, dePar, S=vCov, models=model, type = "general",
               pVal = T, df=Inf)
dd2 <- dd;dd2$x <- as.factor(dd$x)
fit <- glm(y~x-1, family = binomial, data=dd2, weights = n)
mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
summary(mcp)
print(obj, digits = 3)


########################################################################
## unordered values in MCTtest
## placebo-adjusted scale
## two blocks below should give equal results
data(IBScovars)
modlist <- Mods(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1),
                linInt = c(0, 1, 1, 1), doses = c(0, 1, 2, 3, 4))
ancMod <- lm(resp~factor(dose)+gender, data=IBScovars)
drEst <- coef(ancMod)[2:5]
vc <- vcov(ancMod)[2:5, 2:5]
doses <- 1:4
fitMod(doses, drEst, S=vc, model = "sigEmax", placAdj=TRUE, type = "general")
MCTtest(doses, drEst, S = vc, models = modlist, placAdj = TRUE,
        type = "general", df = Inf)


ord <- c(3,4,1,2)
drEst2 <- drEst[ord]
vc2 <- vc[ord,ord]
doses2 <- doses[ord]
fitMod(doses2, drEst2, S=vc2, model = "sigEmax", placAdj=TRUE, type = "general")
MCTtest(doses2, drEst2, S = vc2, models = modlist, placAdj = TRUE,
        type = "general", df = Inf)

## unadjusted scale
## two blocks below should give equal results
ancMod <- lm(resp~factor(dose)-1, data=IBScovars)
drEst <- coef(ancMod)
vc <- vcov(ancMod)
doses <- 0:4
fitMod(doses, drEst, S=vc, model = "sigEmax", type = "general")
MCTtest(doses, drEst, S = vc, models = modlist, type = "general", df = Inf)


ord <- c(3,4,1,2,5)
drEst2 <- drEst[ord]
vc2 <- vc[ord,ord]
doses2 <- doses[ord]
fitMod(doses2, drEst2, S=vc2, model = "sigEmax", type = "general")
MCTtest(doses2, drEst2, S = vc2, models = modlist, type = "general", df = Inf)
