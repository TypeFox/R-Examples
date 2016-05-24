# chapter 4
expon <- function(t,p,x) {
 # p is growth or decay rate in exponential model
 dx <- p[1]*x
 return(dx)
}

expon.rand <- function(t,p,x){
 mu<-p[1]; sd<-p[2] # rate in exponential model
 dx <- rnorm(1,mu,sd)*x
 return(dx)
}

# chapter 6
logistic <- function(t,p,x){
 r<-p[1]; K<-p[2]
 dx <- r*x*(1-x/K)
 return(dx)
}

M3 <- function(t,p,x){
 # p[1] negative for decay
 Kmax<-p[1]; Kh<-p[2]
 dx <- Kmax*x/(x+Kh)
 return(dx)
}

monod.batch <- function(t,p,x){ 
 # x[1] substrate, x[2] microbe
 Kmax<-p[1]; Kh<-p[2]; y<-p[3]; D<-p[4]
 dx1 <- -x[2]*Kmax*x[1]/(x[1]+Kh)
 dx2 <- -dx1*y - D*x[2]
 dx<-c(dx1, dx2)
 return(dx) 
}

# Chapter 7
expon.forced <- function(t,p,x){
 r <- p[1]; Hp <- p[2]; Ha <- p[3]
 dx <- r*x + Hp*x + Ha
 return(dx)
 }

logistic.forced <- function(t,p,x){
 r<-p[1]; K<-p[2]; Hp <- p[3]; Ha <- p[4]
 dx <- r*x*(1-x/K) + Hp*x + Ha
 return(dx)
 }

expon.var <- function(t,p,x){
 # dx = (rm + rd*t + ra*sin(2pi t/T))*x
 r <- p[1] + p[2]*t + p[3]*sin(2*pi*t/p[4])
 dx <- r*x
 return(dx)
 }

logistic.var <- function(t,p,x){
 # K(t) = K0 + Kd*t + Ka*sin(2pi t/T)
 r<- p[1];K <- p[2] + p[3]*t + p[4]*sin(2*pi*t/p[5])
 dx <- r*x*(1-x/K)
 return(dx)
 }

expon.z <- function(t,p,x){
 # periodic impulse train
 tz <- seq(t[1],t[length(t)],p[2])
 return(tz)
}

expon.g <- function(t,p,x,tz){
 # linear disturbance at impulse times
 u <- 0
 for(j in 1:length(tz))
 if(x>0) if(abs(t - tz[j])<0.00001) u <- p[3]*x + p[4]
 return(u)
}

logistic.z <- function(t,p,x){
 # periodic impulse train
 tz <- seq(t[1],t[length(t)],p[3])
 return(tz)
}

logistic.g <- function(t,p,x,tz){
 # linear disturbance at impulse times
 u <- 0
 for(j in 1:length(tz))
 if(x>0) if(abs(t - tz[j])<0.00001) u <- p[4]*x + p[5]
 return(u)
}

#chapter 10
two.stage.cont <- function(t,p,x){ 
 # x[1] juveniles, x[2] mature adults
 # b, m, d1, d2 are p[1:4]
 dx1 <- p[1]*x[2] - p[2]*x[1] - p[3]*x[1]
 dx2 <- p[2]*x[1] - p[4]*x[2]
 dx <- c(dx1,dx2) 
 return(dx)
}

two.stage.cont.delay <- function(t,p,x){ 
 #x[1] juveniles, x[2] mature adults
 #x[3],...,x[n+2] delay vars
 # b, m, d1, d2 are p[1:4]
 n <- p[5]; p[3] <- p[3]/(n+1); y <- array()
 y[1] <- p[1]*x[2] - p[2]*x[1] - p[3]*x[1]
 y[3] <- p[2]*(x[1] - x[3]) - p[3]*x[3]
 for(i in 4:(n+2)) 
 y[i] <- p[2]*(x[i-1] - x[i]) - p[3]*x[i]
 y[2] <- p[2]*x[n+2] - p[4]*x[2]
 dx <- y 
 return(dx)
}

two.stage.x0 <- function(p,X0){
 x0 <- c( X0[1]/(p[5]+1), X0[2], rep(X0[1]/(p[5]+1),p[5]) )
 return(x0)
}

# chapter 11
toxi1.bioacc <- function(t,p,x){
 #x body concentration
 #"ex","ew","V","W","Cw" p[1:5]
 dx <- (p[2]*p[3]/p[4])*p[5] - p[1]*x 
 return(dx)
} 

toxi2.bioacc<- function(t,p,x){ 
 # x[1], x[2], compartmnents
 # 1"ex",2"u",3"Cw",4"k12",5"k21"
 k <- matrix(ncol=2,nrow=2) # transfers
 k[1,1] <- 0; k[1,2] <- p[4]
 k[2,1] <- p[5]; k[2,2] <- 0
 k[1,1] <- -sum(k[,1]); k[2,2] <- -sum(k[,2]) 
 u <- c(p[2], 0.0)*p[3]
 e <- c(0.0,p[1])
 dX <- c(k%*%x + u - e*x)
 return(dX)
 }

toxi.multi.bioacc <- function(t,p,x){ 
 # x[1], x[2], x[3] compartmnents
 # parameters 1"ex",2"u",3"Cw",4"k12",5"k13",6"k21",
 # 7"k23",8"k31",9"k33",10"Kmax",11"Kh",12"Te",13"Du"
 
 # pulse
 if(p[12]>0){
 j <- floor(t/p[12])+1
 if(t< p[12]*(j+p[13]-1)) u1 <- p[3] 
 if(t>=p[12]*(j+p[13]-1)) u1 <- 0
 } else u1 <- p[3]
 u <- c(p[2], 0.0, 0.0)*u1

 k <- matrix(ncol=3,nrow=3) # transfers
 k[1,1] <- 0; k[1,2] <- p[4]; k[1,3] <- p[5]
 k[2,1] <- p[6]; k[2,2] <- 0; k[2,3] <- p[7]
 k[3,1] <- p[8]; k[3,2] <- p[9]; k[3,3] <- 0
 k[1,1] <- -sum(k[,1]); k[2,2] <- -sum(k[,2]); k[3,3] <- -sum(k[,3]) 

 # excretion and breakdown
 e <- c(p[1], 0.0, 0.0)
 b <- c(0.0, p[10]/(p[11]+x[2]), 0.0)
 dX <- c(k%*%x + u -(e+b)*x)
 for(i in 1:3) if(x[i] < 0) x[i] <- 0
 return(dX)
}

# chapter 12
LVint.2sp <- function(t,p,x){ 
 # write model with all positive coefficients
 # let signs of parameters values determine the type of interaction
 r <- p[1:2]; A <- matrix(p[3:6],ncol=2,byrow=T); u <- p[7:8] 
 dX <- c(r +A%*%x)*x +u
 return(dX)
}

holling.tanner <- function(t,p,x){ 
 # x[1] resource, x[2] consumer
 # p[1] r1, p[2] r2, p[3] a11, p[4] a22, p[5] Kmax, p[6] Kh, p[7] h 
 dX1 <- (p[1] - p[3]*x[1]-x[2]*p[5]/(x[1]+p[6]))*x[1] + p[7]
 dX2 <- (p[2] - p[4]*x[2]/x[1])*x[2] + p[7]
 dX<-c(dX1, dX2)
 return(dX)
}

LVint.3sp <- function(t,p,x){ 
 r <- p[1:3]; A <- matrix(p[4:12],ncol=3,byrow=T); u <- p[13:15]
 dX <- (r + A%*%x)*x + u
 return(dX)
}

succession <- function(t,p,x){ 
 k <- matrix(ncol=5,nrow=5) # transfers
 k[1,1:5] <- 0
 k[2,1] <- p[1]; k[2,2:5] <- 0
 k[3,1] <- p[2]; k[3,2] <- p[3]; k[3,3:5] <- 0
 k[4,1:2] <- 0.0; k[4,3] <- p[4]; k[4,4:5] <- 0
 k[5,1:3] <- 0.0; k[5,4] <- p[5]; k[5,5] <- 0
 for(i in 1:5) k[i,i] <- -sum(k[,i]) 
 dX <- c(k%*%x)
 return(dX)
}

succession.z <- function(t,p,x) {
 # periodic impulse train
 tz <- seq(t[1],t[length(t)],p[7])
 return(tz)
}

succession.g <- function(t,p,x,tz){
 u <- rep(0,5)
 for(j in 1:length(tz))
 if(abs(t - tz[j])<0.00001) {
  u[2:5] <- rep(p[6],4)*x[2:5]
  u[1] <- -sum(u[2:5])
 }
 return(u)
}

# chapter 13

nut.cycle <- function(t,p,x){ 
 k <- matrix(ncol=4,nrow=4) # transfers
 k[,] <- 0
 k[1,3] <- p[1]; k[1,4] <- p[2]; k[2,1] <- p[3]
 k[3,2] <- p[4]; k[4,3] <- p[5]
 for(i in 1:4) k[i,i] <- -sum(k[,i]) 
 U <- c(sum(p[8:9]),0,0,0) # abs input
 uc <- c(0,0,0,p[7]) # prop input
 dX <- c(k%*%x + uc*x +U)
 return(dX)
}

nut.cycle.z <- function(t,p,x){
 # periodic impulse train
 tz <- seq(t[1],t[length(t)],p[11])
 return(tz)
}

nut.cycle.g <- function(t,p,x,tz){
 u <- rep(0,4)
 for(j in 1:length(tz))
 if(abs(t - tz[j])<0.00001){
 u[2] <- p[6]*x[2]
 u[1] <- p[10] 
 }
 return(u)
}

#chapter 14
DO.PP.pond <- function(t,p,x){
 # parameters
 D <- p[1]; T1 <- p[2]; dT <- p[3]
 k <- p[4]; zd <- p[5]; Pmax <- p[6]
 alpha <- p[7]; lat <- p[8]; day1 <- p[9]
 Lm <- p[10]; sdr <- p[11]; Rsp <- p[12]
 # determine day
 i <- floor(t/24); day <- day1+i
 hr <- t -24*i -12 # hr within the day
 temp <- T1 + dT*t # temp linear trend
 Xs <- 14.6 - ((14.6-8.6)/25)*temp # DO saturation
 Ls <- sun.rad.hr(day,alpha,hr,Lm,sdr)
 dX <- D*(Xs-x) + PPT.Smith(Ls,k,zd,Pmax,alpha) - Rsp
 return(dX)
}

nut.river <- function(t,p,x){
 dX <- x # define array
 V <- p[1]; Q<- p[2]; W <- p[3] # vol, flow, discharge
 Pmax <- p[4]; Kh <- p[5] # uptake M3 parameters
 X1.u <- p[6]; X2.u <- p[7] # boundary conditions
 # x[1] nitrogen x[2] algae
 W.V <- (W/V)*(24*60*60) # loading convert s to d
 Q.V <- (Q/V)*(24*60*60) # depuration convert s to d
 U2 <- Pmax*(x[1]/(x[1]+Kh))*x[2] # M3 uptake
 dX[1] <- W.V + Q.V*(X1.u - x[1]) - U2 
 dX[2] <- Q.V*(X2.u - x[2]) + U2
 return(dX)
}

# chapter 15

soilwat.1 <- function(t,p,x){
 # infilt rate capacity linear with deficit
 # infiltration parameters in mm/hr 
 fd <- p[1]; Ks <- p[2]
 # soil capacity in mm, field capacity, porosity
 Z <- p[3]; Fc.coeff <- p[4]; porosity <- p[5]
 rain.int <- p[6]; rain.dur <- p[7] 
 # soil sat capacity and field capacity
 cap.soil <- Z * porosity
 Fc <- cap.soil * Fc.coeff
 # using deficit to determine infiltration cap
 deficit <- cap.soil - x
 # slope or linear component mm/hr per mm
 linear <- (fd-Ks)/cap.soil
 # infiltration cap in mm/hr
 infilt <- Ks + deficit*linear
 if (infilt > fd) infilt <- fd
 if (x >= cap.soil) f<- Ks
 else f<- infilt
 # compare to rain 
 if(t <= rain.dur) rain <- rain.int else rain <- 0
 if(rain< f ) q <- rain else q <- f
 # percolation
 if(x >= Fc) g <- Ks * (x-Fc)/(cap.soil-Fc)
 else g <- 0
 dx <- q-g
 return(dx)
}

green.ampt.maila <- function(t,p,x){
 Pf<- p[1]; Ks<- p[2] 
 porosity <- p[3]; init <- p[4]; dtr <- p[5]
 if(x <= 0) x <- 10^-5
 deficit <- porosity - init
 num <- 2*dtr*Ks*(1+ Pf*deficit/x)
 den <- 2 - dtr*(-Ks*Pf*deficit/x^2)
 dx <- num/den
 return(dx)
}

green.ampt.ramos <- function(t,p,x){
 Pf<- p[1]; Ks<- p[2] 
 porosity <- p[3]; init <- p[4]; Ic <- p[5]
 if(x <= 0) x <- 10^-5
 deficit <- porosity - init
 dx<- Ks*(1+ Pf*deficit/x)
 df <- -Ks*Pf*deficit/x^2
 return(c(dx,df))
}

sim.ga <- function(t,p){
 x <- array(); f <- x; q <- x; Q <- x 
 x[1] <- 10^-5; f[1] <- 0; q[1] <- 0; Q[1]<-0
 dtr <- t[2]-t[1] 
 for(i in 2:length(t)){
 # integration
 dF <- green.ampt.maila(t[i],p,x[i-1])
 x[i] <- x[i-1] + dF 
 # infiltration capacity
 f[i] <- dF/dtr
 rain <- p[6] 
 # actual rate and infiltrated depth
 if(f[i] < rain) q[i] <- f[i] else q[i] <- rain
 Q[i] <- Q[i-1]+q[i]*dtr
 } 
 mat <- matrix(1:2, c(2,1),byrow=T)
 layout(mat,rep(7,2),rep(7/2,2),respect=TRUE)
 par(mar=c(4,4,1,.5), xaxs="r", yaxs="r")

 matplot(t,cbind(x,Q),type="l",col=1,ylab="Depth [mm]",xlab="Time [h]",ylim=c(0,50))
 legend("topleft",lty=1:2,col=1,legend=c("Capacity", "Actual"))

 matplot(t,cbind(f,q),type="l",col=1,ylab="Rate [mm/h]",xlab="Time [h]",ylim=c(0,500))
 legend("topright",lty=1:2,col=1,legend=c("Capacity", "Actual"))

 out <- data.frame(t,x,f,q,Q)
 return(round(out,2))
} # end of sim.ga



soilwat.n <- function(t,p,x){
 # infilt rate capacity linear with deficit
 # infiltration parameters in mm/hr 
 fd <- p[1]; Ks <- p[2]
 # soil capacity in mm, field capacity, porosity
 Z <- p[3]; Fc.coeff <- p[4]; porosity <- p[5]
 rain.int <- p[6]; rain.dur <- p[7]
 n <- length(x)
 # soil sat capacity and field capacity
 cap.soil <- Z * porosity
 Fc <- cap.soil * Fc.coeff
 # using deficit to determine infiltration cap
 deficit <- cap.soil - x[1]
 # slope or linear component mm/hr per mm
 linear <- (fd-Ks)/cap.soil
 # infiltration cap in mm/hr
 infilt <- Ks + deficit*linear

 if (infilt > fd) infilt <- fd
 if (x[1] >= cap.soil) f<- Ks
 else f<- infilt

 if(t <= rain.dur) rain <- rain.int else rain <- 0
 if(rain< f ) q <- rain else q <- f
 g <- array(); dx <- array()

 # percolation from each layer (except last)
 for(i in 1:(n-1)){ 
 if(x[i] >= Fc) {
  if(x[i+1]>=cap.soil) g[i] <- Ks 
  else g[i] <- Ks * (x[i]-Fc)/(cap.soil-Fc)
 } else g[i] <- 0
 if(i==1) dx[i] <- q-g[i] else dx[i] <- g[i-1]-g[i]
 } # end of for

 # last layer
 if(x[n] >= Fc) g[n]<- Ks * (x[n]-Fc)/(cap.soil-Fc)
 else g[n] <- 0
 dx[n] <- g[n-1]-g[n]
 
 return(dx)
}



