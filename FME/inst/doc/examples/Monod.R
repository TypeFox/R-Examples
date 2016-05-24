##==============================================================================
## Fitting a Monod function to data
## example as in Berthoux and Brown, 2002.
##==============================================================================


library(FME)

# 1. the observations
#---------------------
Obs <- data.frame(x=c(   28,  55,   83,  110,  138,  225,  375),   # mg COD/l
                  y=c(0.053,0.06,0.112,0.105,0.099,0.122,0.125))   # 1/hour
Obs2 <- cbind(Obs, err= Obs$y*0.1)

plot(Obs,pch=16,cex=2,xlim=c(0,400),ylim=c(0,0.15),
     xlab="mg COD/l",ylab="1/hr")

# 2. the Monod model
#---------------------
Model <- function(p,x)   return(data.frame(x=x,y=p[1]*x/(x+p[2])))

# 3. Fitting the model to the data
#---------------------
# define the residual function
Residuals  <- function(p) (Obs$y-Model(p,Obs$x)$y)  #... model residuals

# fit the model to data
P      <- modFit(f=Residuals,p=c(0.1,1))

# plot best-fit model
x      <-0:375
lines(Model(P$par,x))

# Now with model error...
mcost <- function(p) {
  Mod <- Model(p,Obs$x)
  modCost(Mod,Obs2,x="x",err="err")
}

# fit the model to data
P2     <- modFit(f=mcost,p=P$par)
lines(Model(P2$par,x),col="red")


# summary of fit
sP    <- summary(P)
sP[]
print(sP)

# 4. MCMC analysis
#---------------------
# estimate of parameter covariances (to update parameters) and the model variance
Covar   <- sP$cov.scaled * 2.4^2/2
s2prior <- sP$modVariance

# set wvar0 = 0.2 to avoid too much updating model variance
MCMC <- modMCMC(f=Residuals,p=P$par,jump=Covar,niter=3000,
                var0=s2prior,wvar0=1,updatecov=100,lower=c(0,0))

plot(MCMC,Full=TRUE)
pairs(MCMC)
cor(MCMC$pars)
cov(MCMC$pars)
sP$cov.scaled

sR<-sensRange(parInput=MCMC$pars,func=Model,x=1:375)
plot(summary(sR))
points(Obs)
hist(MCMC,Full=TRUE,col="darkblue")

# increase acceptance rate with delayed rejection
MCMC2 <- modMCMC(f=Residuals,p=P$par,jump=Covar,niter=3000, ntrydr=3,
                var0=s2prior,wvar0=1,updatecov=100,lower=c(0,0))
plot(MCMC2,Full=TRUE)
pairs(MCMC2)
MCMC2$count
# show only the posterior model error and the function value
hist (MCMC,Full=TRUE,which=NULL,col="darkgreen")

