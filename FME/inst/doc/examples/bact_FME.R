## =============================================================================
## A simple model of bacteria, growing on a substrate
## in a batch culture
## Model from Soetaert and Herman, 2009
## Checked 24-08-09
## =============================================================================

mf <- par(mfrow=c(2,2))
require(FME)

##------------------------------------------------------------------------------
##                           the Model
##------------------------------------------------------------------------------

pars <- list(gmax =0.5,eff = 0.5,
              ks =0.5, rB =0.01, dB =0.01)

solveBact <- function(pars, tout = seq(0,50,by=0.5)) {
  derivs <- function(t,state,pars) {    # returns rate of change
    with (as.list(c(state,pars)), {
      dBact = gmax*eff*Sub/(Sub+ks)*Bact - dB*Bact - rB*Bact
      dSub  =-gmax    *Sub/(Sub+ks)*Bact + dB*Bact
      return(list(c(dBact,dSub)))
    })
  }
  state   <- c(Bact=0.1,Sub = 100)
  ## ode solves the model by integration...
  return(as.data.frame(ode(y=state,times=tout,func=derivs,parms=pars)))
}

out <- solveBact(pars)

plot(out$time,out$Bact,ylim=range(c(out$Bact,out$Sub)),
     xlab="time, hour",ylab="molC/m3",type="l",lwd=2)
lines(out$time,out$Sub,lty=2,lwd=2)
lines(out$time,out$Sub+out$Bact)

legend("topright",c("Bacteria","Glucose","TOC"),
       lty=c(1,2,1),lwd=c(2,2,1))

##------------------------------------------------------------------------------
##       Local sensitivity analysis : sensitivity functions
##------------------------------------------------------------------------------

## sensitivity functions
SnsBact <- sensFun(func=solveBact,parms=pars,
                   sensvar="Bact",varscale=1)
head(SnsBact)
plot(SnsBact)
plot(SnsBact,which="Bact")

summary(SnsBact)
plot(summary(SnsBact))

SF<- sensFun(func=solveBact,parms=pars,
             sensvar=c("Bact","Sub"),varscale=1)
head(SF)
tail(SF)
plot(SF,legpos="bottomright")

pairs(SF, which=c("Bact","Sub"))
pairs(SF, which=c("Sub","Bact"))

summary(SF,var=TRUE)

## Bivariate sensitivity
##-----------------------------
pairs(SnsBact)
mtext(outer=TRUE,side=3,line=-2,
      "Sensitivity functions",cex=1.5)

## pairwise correlation
cor(SnsBact[,-(1:2)])

## multivariate sensitivity analysis
##-----------------------------
Coll <- collin(SnsBact[,-(1:2)])

## The larger the collinearity, the less identifiable the data set
Coll

nc <- ncol(Coll)
plot(Coll,log="y")

## 20 = magical number above which there are identifiability problems
abline(h=20,col="red")
Coll [Coll[,"collinearity"]<20&Coll[,"N"]==4,]


##------------------------------------------------------------------------------
##       Global sensitivity analysis : Sensitivity ranges
##------------------------------------------------------------------------------

## the sensitivity parameters
parRanges <- data.frame(min=c(0.4,0.4,0.),max=c(0.6,0.6,0.02))
rownames(parRanges)<- c("gmax","eff","rB")
parRanges

tout    <- 0:50
## sensitivity to rB; equally-spaced parameters ("grid")
SensR <- sensRange(func=solveBact,parms=pars,dist="grid",
                   sensvar="Bact",parRange=parRanges[3,],num=50)

Sens  <-summary(SensR)
plot(Sens,legpos="topleft",xlab="time, hour",ylab="molC/m3",
     main="Sensitivity to rB")

## sensitivity to all; latin hypercube
Sens2 <- summary(sensRange(func=solveBact,parms=pars,dist="latin",
                   sensvar=c("Bact","Sub"),parRange=parRanges,num=50))
plot(Sens2,xlab="time, hour",ylab="molC/m3",
     main="Sensitivity to gmax,eff,rB")

##------------------------------------------------------------------------------
##                Fitting the model to the data
##------------------------------------------------------------------------------

## the "data"
Data <- matrix (nc=2,byrow=TRUE,data=
c(  2,  0.14,    4,  0.2,     6,  0.38,    8,  0.42,
   10,  0.6,    12,  0.107,  14,  1.3,    16,  2.0,
   18,  3.0,    20,  4.5,    22,  6.15,   24,  11,
   26, 13.8,    28, 20.0,    30,  31 ,    35, 65, 40, 61)
)
colnames(Data) <- c("time","Bact")
head(Data)

Data2 <- matrix(nc=2,byrow =TRUE,data=c(2,100,20,93,30,55,50,0))
colnames(Data2) <- c("time","Sub")

## Objective function to minimise; parameters gmax and eff are fitted

Run <- function(x) {
 pars[c("gmax","eff")]<- x
 solveBact(pars)
}


Objective <- function (x,out=Run(x)) {   # Model cost
  Cost <- modCost(obs=Data2,model=out)   # observed data in 2 data.frames
  return(modCost(obs=Data,model=out,cost=Cost))
}

sF <- sensFun(func=Objective,parms=pars[c("gmax","eff")],varscale=1)
collin(sF)

## 2. modFit finds the minimum; parameters constrained to be > 0
print(system.time(Fit<-modFit(p=c(0.5,0.5),f=Objective)))

Fit
summary(Fit)

## Run best-fit model....
out <- Run(Fit$par)

## Model cost
Cost <- Objective(Fit$par,out)    # observed data in 2 data.frames

## Plot residuals
plot(Cost,xlab="time, hour",ylab="molC/m3",main="residuals",cex=1.5)


## plot output
plot(out$time,out$Bact,ylim=range(out$Bact),
     xlab="time, hour",ylab="molC/m3",type="l",lwd=2)
points(Data,cex=2,pch=18)

##------------------------------------------------------------------------------
##                        MCMC application
##------------------------------------------------------------------------------

## estimate of parameter covariances (to update parameters) and the model variance
sP <- summary(Fit)
Covar   <- sP$cov.scaled * 2.4^2/2
s2prior <- sP$modVariance

## set nprior = 2 to avoid too much updating model variance
MCMC <- modMCMC(f=Objective,p=Fit$par,jump=Covar,niter=1000,
                var0=s2prior,wvar0=1,updatecov=10)

plot(MCMC,Full=TRUE)
pairs(MCMC)
cor(MCMC$pars)
cov(MCMC$pars)
sP$cov.scaled

sR<-sensRange(parInput=MCMC$pars,func=Run)
plot(summary(sR))
points(Data2)

## Use delayed rejection
MCMC2 <- modMCMC(f=Objective,p=Fit$par,jump=Covar,niter=1000,
                var0=s2prior,wvar0=1,updatecov=10, ntrydr=3)
pairs(MCMC2)
cor(MCMC2$pars)


## use functions from the coda package
MC <-as.mcmc(MCMC$pars)
MC2<-as.mcmc(MCMC2$pars)

MClist <- as.mcmc.list(list(MC,MC2))
gelman.diag(MClist)

cumuplot(MC2)
