## ----eval=FALSE----------------------------------------------------------
#  G1=function(t){-sin(2*pi*t)}
#  Q0=function(t){1}

## ----fig.align = 'center'------------------------------------------------
Xs <- 0                 # Initial state
Xt <- seq(-3/2,3/2,1/50)# Possible future states
s  <- 0                # Starting time
t  <- 1                # Final time
mu    <- 0.5             # Drift parameter
sigma <- 0.25            # Diffusion coefficient

# True transition density
plot(dnorm(Xt,Xs+mu*(t-s),sigma*sqrt(t-s))~Xt,main='Transition density',type='l')

## ----fig.align = 'center'------------------------------------------------
library(DiffusionRgqd) 

# Remove any existing coefficients:
GQD.remove()    

# Define the model coefficients:
G0 <- function(t){mu}
Q0 <- function(t){sigma^2}

# Calculate the transitional density:
BM <- GQD.density(Xs,Xt,s,t)

# Plot the transitional density:
plot(dnorm(Xt,Xs+mu*(t-s),sigma*sqrt(t-s))~Xt,main='Transition density',type='l')
lines(BM$density[,100]~BM$Xt,col='blue',lty='dashed',lwd=2)

## ------------------------------------------------------------------------
 
# Remove any existing coefficients
GQD.remove()         

a = 0.5; b =5; sigma =0.35; # Some parameter values.

# Define drift Coefficients.
G0 <- function(t){a*b}    
G1 <- function(t){-a}
# Define sinusoidal diffusion coefficient.
Q1 <- function(t){sigma^2}

states     <-  seq(1,9,1/10)  # State values
initial    <-  6              # Starting value of the process
Tmax       <-  5              # Time horizon
Tstart     <-  1              # Time starts at 1
increment  <-  1/100          # Incremental time steps

# Generate the transitional density
M <- GQD.density(Xs=initial,Xt=states,s=Tstart,t=Tmax,delt=increment)


## ----fig.align = 'center'------------------------------------------------
persp(x=M$Xt,y=M$time,z=M$density,col=3,xlab='State (X_t)',ylab='Time (t)',
      zlab='Density f(X_t|X_s)',border=NA,shade=0.5,theta=145)

## ----fig.align = 'center'------------------------------------------------
# What is returned:
names(M)

# Plot the mean trajectory:
plot(M$cumulants[1,]~M$time,col='blue',main='Mean trajectory',type='l')

## ----fig.align = 'center'------------------------------------------------
# Plot the variance trajectory:
plot(M$cumulants[2,]~M$time,col='blue',main='Variance trajectory',type='l')

## ----eval=FALSE----------------------------------------------------------
#  browseVignettes('DiffusionRgqd')

## ----eval=FALSE----------------------------------------------------------
#  demo(package = 'DiffusionRgqd')

