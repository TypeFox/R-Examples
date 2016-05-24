##################################################################################
###                                                                            ###
###  Example script for First Passage Time simulation for Gaussian diffusions  ###
###                                                                            ###
###       DIFFUSION WITH TIME-DEPENDENT COEFFS THROUGH A LINEAR BOUNDARY       ###
###                                                                            ###
##################################################################################

library(GaDiFPT)


#############################################################
######         FIRST PART: BUILDING THE  PROCESS       ######
#############################################################

cat('#################################################################### \n')
cat('#####             First Passage Time Simulation for            ##### \n')
cat('#####    a diffusion process with time-dependent coefficients  ##### \n')
cat('#####                 through a constant boundary              ##### \n')
cat('#################################################################### \n \n \n')


Logistic <- diffusion(c("-x/theta(t) + rho(t)/theta(t) + mu","sigma2"))

print(Logistic)


############################################################
####                   INPUT AND SETUP                  ####
############################################################


# Input the infinitesimal moments of the process
# along with the related parameters.


erre <- 4.0
kappa <- 1.0
theta <- 40 
rho <- 0.0
mu <- 0.25 
sigma2 <- 0.2 


# Input the threshold parameters

Scost <- 15.5
Sslope <- 0.0
Stype <- "constant"


# Input the user-provided parameters (initial time and state,
# final time, initial timestep...).  


# Initial time t0 (in milliseconds): 
t0 <- 0.0

# Rest potential x0: 
x0 <- 7.5
#x0 <- 0

# Final time Tfin (in milliseconds): 
Tfin <- 4000 

# Initial timestep delta t (in milliseconds): 
deltat <- 0.5


# tolerance for the stopping criterion in the refinement procedure
tol1 <- 0.001

# number of timesteps
N <- floor((Tfin - t0)/deltat)  


# number of spikes to be simulated

M <- 10000

# flag for requiring numerical integration 
# quadflag = 1 for integration up to Tfin
# quadflag = 0 for integration just to the time when
# the stationary value of the hazard function is reached 

quadflag <- 1

# temporary flag for execution in RStudio
RStudioflag <- TRUE


param <- inputlist(mu,sigma2,Stype,t0,x0,Tfin,deltat,M,quadflag,RStudioflag)

print(param)

## Specify a filename for writing the output
fileout <- "results_Logistic.out"


# building the infinitesimal moments

aaa <- function(t) {
  temp <- erre*(t-t0)/theta
  aaa <- -(1.0 + kappa*exp(-temp))/theta
}

bbb <- function(t) {
  temp <- erre*(t-t0)/theta
  bbb <- rho*(1.0 + kappa*exp(-temp))/theta
  bbb <- bbb + mu
}


# building the threshold as S(t) and threshold derivative Sp(t). 

SSS <- function(t) {
  SSS <- Scost + Sslope*t
}

SSSp <- function(t) {
  SSSp <- Sslope 
}

####       INITIALIZATION OF VECTORS

# vector of timesteps
tempi <- numeric(N+1)

# mean of the process
mp <- numeric(N+1)

# covariance of the process: cov = up*vp
up <- numeric(N+1)
vp <- numeric(N+1)

# dummy vector
app <- numeric(N)


####    EVALUATION OF MEAN AND COVARIANCE

tempi <- seq(t0, by=deltat, length=N+1)

# call to the function evaluating the mean mp and the two covariance factors up, vp

dum <- vectorsetup(param)
mp <- dum[,1]
up <- dum[,2]
vp <- dum[,3]


####    PRELIMINARY PLOTS

# produces plots of the threshold and of the mean of the process 
# to check the subthreshold regimen hypothesis
# and gives a preliminary estimate of the time when
# the stationary values are reached

## plot of S and m

splot <- S(tempi)
mp1 <- mp - sqrt(2*sigma2)
mp2 <- mp + sqrt(2*sigma2)
matplot(tempi, cbind(mp,mp1,mp2,splot),type="l",lty=c(1,2,2,1),lwd=1,
        main="mean of the process vs. threshold",xlab="time(ms)",ylab="")
legend("bottomright",c("mean","threshold"),
       lty=c(1,1),col=c("black","blue"))

# first estimate for tmax
Nmax <- which.min(abs(mp[2:(N+1)]-mp[1:N]))


cat('********** \n Estimate FPT density and generate spikes \n') 


  ############################################################
  ####       SECOND PART: EVALUATION OF FPT DENSITY       ####
  ############################################################
  
  # evaluation of g0 and gg0 by numerical integration of Volterra integral equation 
  # up to Tfin
# if a good estimate of the asymptotic firing rate is available then
# density evaluation up to Tfin can be avoided

N1 <- N  

if (quadflag == 0)   N1 <- max(c(Nmax,N/4))
N1p1 <- N1+1


answer <- FPTdensity_byint(param,N1)

plot(answer)


############################################################
####          THIRD PART: SIMULATION OF SPIKES          ####
############################################################


if (M > 0) {
  spikes <- FPTsimul(answer,M)
  
  histplot(spikes,answer)
  
}
############################################################
####          FOURTH PART: SUMMARY OF RESULTS           ####
############################################################


res_summary(answer,M,fileout) 