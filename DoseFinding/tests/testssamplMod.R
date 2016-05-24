require("DoseFinding")
########################################################################
## test Bayesian fitting
## compare bFitMod on example data set with jags

data(biom)

## (i) fit biom data
## data to fit
model <- "sigEmax"
anMod <- lm(resp~factor(dose)-1, data=biom)
drFit <- coef(anMod);y <- drFit
vCov <- vcov(anMod)
Omega <- solve(vCov)#+diag(5)*1000
dose <- sort(unique(biom$dose))
nD <- length(dose)
prior <- list(t = c(0, sqrt(2), 3), norm = c(0, sqrt(2)),
              beta=c(0,1,1,1), lnorm=c(0, sqrt(0.5)))
res <- bFitMod(dose, drFit, vCov, model = "sigEmax", prior=prior, nSim = 100)

## ## jags code (commented out, only for "manual" testing)
## library(rjags)
## path <- "~/Projekte/DoseFindingPackage/"
## modelstr <- "
## model{
##   y[] ~ dmnorm(mu[], Omega[,])
##   for(i in 1:nD){
##     mu[i] <- E0 + (Emax*dose[i]^h)/(dose[i]^h+ED50^h)
##   }
##   E0 ~ dt(0, 0.5, 3)
##   Emax ~ dnorm(0, 0.5)
##   ED50 ~ dunif(0,1)
##   h ~ dlnorm(0, 2)
## }
## "
## file <- paste(path, "mod.txt", sep="")
## cat(modelstr, file = file)

## ## data
## jags.data <- list(y=y, nD=nD, dose=dose, Omega=Omega)
## jags.inits <-  list("E0"=0,"Emax"=0,"ED50"=0.5,"h"=1)
## mod <- jags.model(file, jags.data, jags.inits, n.chains = 3)
## samp <- jags.samples(mod, c("E0","Emax","ED50", "h"), 100000)

## quantile(samp$E0, c(0.05,0.25,0.5,0.75,0.95))
## quantile(res$samples[,1], c(0.05,0.25,0.5,0.75,0.95))
## quantile(samp$Emax, c(0.05,0.25,0.5,0.75,0.95))
## quantile(res$samples[,2], c(0.05,0.25,0.5,0.75,0.95))
## quantile(samp$ED50, c(0.05,0.25,0.5,0.75,0.95))
## quantile(res$samples[,3], c(0.05,0.25,0.5,0.75,0.95))
## quantile(samp$h, c(0.05,0.25,0.5,0.75,0.95))
## quantile(res$samples[,4], c(0.05,0.25,0.5,0.75,0.95))

## cor(cbind(as.numeric(samp$E0[,,]),
##           as.numeric(samp$Emax[,,]),
##           as.numeric(samp$ED50[,,]),
##           as.numeric(samp$h[,,])))
## cor(res$samples)

## (ii) now run with inflated variance (essentially sample prior)
vCov <- vcov(anMod)*100000
Omega <- solve(vCov)#+diag(5)*1000
res <- bFitMod(dose, drFit, vCov, model = "sigEmax", prior=prior, nSim = 100)

## ## jags code (commented out, only for "manual" testing)
## jags.data <- list(y=y, nD=nD, dose=dose, Omega=Omega)
## mod <- jags.model(file, jags.data, jags.inits, n.chains = 3)
## samp <- jags.samples(mod, c("E0","Emax","ED50", "h"), 100000)

## quantile(samp$E0, c(0.05,0.25,0.5,0.75,0.95))
## quantile(res$samples[,1], c(0.05,0.25,0.5,0.75,0.95))
## quantile(samp$Emax, c(0.05,0.25,0.5,0.75,0.95))
## quantile(res$samples[,2], c(0.05,0.25,0.5,0.75,0.95))
## quantile(samp$ED50, c(0.05,0.25,0.5,0.75,0.95))
## quantile(res$samples[,3], c(0.05,0.25,0.5,0.75,0.95))
## quantile(samp$h, c(0.05,0.25,0.5,0.75,0.95))
## quantile(res$samples[,4], c(0.05,0.25,0.5,0.75,0.95))

## cor(cbind(as.numeric(samp$E0[,,]),
##           as.numeric(samp$Emax[,,]),
##           as.numeric(samp$ED50[,,]),
##           as.numeric(samp$h[,,])))
## cor(res$samples)

########################################################################
## test bootstrap fitting

vCov <- vcov(anMod)
bnds <- matrix(c(0.001, 0.5, 1.5, 10), 2, 2)
res <- bFitMod(dose, drFit, vCov, model = "sigEmax", nSim = 100, bnds=bnds,
               type = "bootstrap")
dd <- dose[-1];resp <- drFit[2:5]-drFit[1]
vc <- cbind(-1,diag(4))%*%vCov%*%t(cbind(-1,diag(4)))
res <- bFitMod(dd, resp, vc, model = "linear", nSim = 100, bnds=bnds,
               placAdj = TRUE, type = "bootstrap")

########################################################################
## test dose calculations, when model = "linInt" and placAdj=TRUE

data(IBScovars)
anovaMod <- lm(resp~factor(dose)+gender, data=IBScovars)
drFit <- coef(anovaMod)[2:5] # placebo adjusted estimates at doses
vCov <- vcov(anovaMod)[2:5,2:5]
dose <- sort(unique(IBScovars$dose))[-1]
fm <- fitMod(dose, drFit, S=vCov, model = "linInt", type = "general", placAdj=TRUE)
ED(fm, 0.25)
ED(fm, 0.5)
ED(fm, 0.75)
ED(fm, 0.95)
TD(fm, 0.2)
TD(fm, 0.3)
TD(fm, 0.4)

prior <- list(norm = c(0,1000), norm = c(0,1000),
              norm = c(0,1000), norm = c(0,1000))
gsample <- bFitMod(dose, drFit, vCov, model = "linInt", placAdj=TRUE,
                   start = c(1, 1, 1, 1), nSim = 1000, prior = prior)
td1 <- TD(gsample, 0.3)
td2 <- TD(gsample, 0.3, TDtype="d", doses = seq(0,4,length=101))
ed1 <- ED(gsample, 0.8)
ed2 <- ED(gsample, 0.8, EDtype="d", doses = seq(0,4,length=101))

