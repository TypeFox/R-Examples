##################################################################################
###                                                                            ###
###  Example script for First Passage Time simulation for Gaussian diffusions  ###
###                                                                            ###
###                 ORNSTEIN-UHLENBECK PROCESS WITH CURRENT                    ###
###                        THROUGH CONSTANT BOUNDARY                           ###                        
##################################################################################


   library(GaDiFPT)


#############################################################
######         FIRST PART: BUILDING THE  PROCESS       ######
#############################################################

cat('#################################################################### \n')
cat('#####               First Passage Time Simulation              ##### \n')
cat('#####            for the Ornstein-Uhlenbeck process            ##### \n')
cat('#####                 with additional current                  ##### \n')
cat('#####               through a constant boundary                ##### \n')
cat('#################################################################### \n \n \n')


OrnUhlCur <- diffusion(c("-x/theta + mu + i0*exp(-(t-t0)/theta1)","sigma2"))

print(OrnUhlCur)


############################################################
####                   INPUT AND SETUP                  ####
############################################################


# Input the infinitesimal moments of the process
# along with the related parameters.


theta <- 20
mu <- 0.1
sigma2 <- 0.25
i0 <- 0.5
theta1 <- 10

# Input the threshold parameters

Scost <- 6
Sslope <- 0
Stype <- "constant"

# Input the user-provided parameters (initial time and state,
# final time, initial timestep...).  


# Initial time t0 (in milliseconds): 
t0 <- 0.0

# Rest potential x0: 
x0 <- 1.0

# Final time Tfin (in milliseconds): 
Tfin <- 500    

# Initial timestep delta t (in milliseconds): 
deltat <- 0.1

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
fileout <- "results_OrnUhlCur.out"


## building the infinitesimal moments
aaa <- function(t) {
  aaa <- -1.0/theta + 0.0*t
}

b <- function(t) {
  bbb <- mu + i0*exp(-(t-t0)/theta1)
}

# building the threshold as S(t) and threshold derivative Sp(t)

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
