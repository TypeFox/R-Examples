## ------------------------------------------------------------------------
 library(DiffusionRgqd)
 GQD.remove()

## ------------------------------------------------------------------------
 G0 <- function(t){2*(10+sin(2*pi*(t-0.5)))}
 G1 <- function(t){-2}
 Q1 <- function(t){0.25*(1+0.75*(sin(4*pi*t)))}

## ------------------------------------------------------------------------
 states    <- seq(5,15,1/10)
 initial   <- 8
 Tmax      <- 5
 Tstart    <- 1
 increment <- 1/100

 M <- GQD.density(Xs=initial,Xt=states,s=Tstart,t=Tmax,delt=increment)

## ----fig.align = 'center'------------------------------------------------
persp(x=M$Xt,y=M$time,z=M$density,col=3,xlab='State (X_t)',ylab='Time (t)',
      zlab='Density f(X_t|X_s)',border=NA,shade=0.5,theta=145)

## ------------------------------------------------------------------------
# Set some parameter values:
a <- 0.5; b <- 0.6; cc <- 1; X0 <- 0.5;

# Give the coefficient form of the diffusion.:
G0 <- function(t){a*(b*(1+0.25*sin(pi*t)))}
G1 <- function(t){-a}
Q1 <- function(t){cc}
Q2 <- function(t){-cc}
# What states should the density be evaluated at:
states=seq(1/1000,1-1/1000,1/1000)

# Generate the transitional density for truncation orders 4,6 and 8:
res1 <- GQD.density(X0,states,0,2,1/100,Dtype='Beta',Trunc=c(4,4))
res2 <- GQD.density(X0,states,0,2,1/100,Dtype='Beta',Trunc=c(6,6),print.output = FALSE)
res3 <- GQD.density(X0,states,0,2,1/100,Dtype='Beta',Trunc=c(8,8),print.output = FALSE)

## ----fig.align='center'--------------------------------------------------
N <- 100000      # Number of trajectories to simulate.
smdelt <- 1/1000 # Simulation stepsize.
d <- 0           # A variable for tracking time.
X <- rep(X0,N)
# Take snapshots of the density at these times:
when <- c(0,0.25,0.5,1,1.75)
for(i in 2:(2/smdelt))
{
# A Euler-Maruyama type approximation:
X <- pmax(pmin(X + (G0(d)+G1(d)*X)*smdelt
+sqrt(Q1(d)*X+Q2(d)*X^2)*rnorm(length(X),sd=sqrt(smdelt)),1),0)
d <- d+smdelt
if(any(when==round(d,3)))
{
index <- which(res1$time==round(d,3))
hist(X,col='#F7F7F7',freq=F,breaks=30,
main=paste0('Transitional Density at t = ',round(d,3)),ylim=c(0,3))
lines(res1$density[,index]~res1$Xt,col='#1B7837',lty='dotdash',lwd=1)
lines(res2$density[,index]~res2$Xt,col='#D92120',lty='solid',lwd=1)
lines(res3$density[,index]~res3$Xt,col='#5289C7',lty='dashed',lwd=2)
legend('top', legend=c('Truncation: (4,4).','Truncation: (6,6).',
'Truncation: (8,8).'), col = c('#1B7837','#D92120','#5289C7'),
lty=c('dotdash','solid','dashed'),lwd=c(1,1,2),bg='white')
}
}


## ----fig.align = 'center'------------------------------------------------
persp(x=res2$Xt,y=res2$time,z=res2$density,col='lightblue',xlab='State (X_t)',ylab='Time (t)', zlab='Density f(X_t|X_s)',border=NA,shade=0.5,theta=145)

## ------------------------------------------------------------------------
library(DiffusionRgqd)
# Remove any existing coefficients:
GQD.remove()
# Define the X dimesnion coefficients:
a10 <- function(t){1.5}
a11 <- function(t){-0.4}
c10 <- function(t){0.05}
# Define the Y dimension coefficients:
b01 <- function(t){-1.5}
b11 <- function(t){0.4}
b02 <- function(t){-0.2}
f01 <- function(t){0.1}
# Approximate the transition density
res <- BiGQD.density(Xs=5,Ys=5,Xt=seq(3,8,length=50),Yt=seq(2,6,length=50),s=0,
t=10,delt=1/100)

## ----fig.align = 'center', warn =FALSE-----------------------------------
# Load simulated trajectory of the joint expectation:
data(SDEsim3)
attach(SDEsim3)
# Record graphs at time points along the trajectory:
time.index <- c(10,200,750,1000) +1
# Make some colour palettes
library(colorspace)
colpal=function(n){rev(sequential_hcl(n,power=1,l=c(40,100)))}
for(i in time.index)
{
# Now illustrate the density using a contour plot:
filled.contour(res$Xt,res$Yt,res$density[,,i],
main=paste0('Transition Density \n (t = ',res$time[i],')'),
color.palette =colpal,
xlab='Prey',ylab='Preditor',plot.axes=
{
# Add trajectory of simulated expectation:
lines(my~mx,col='black',lty='dashed',lwd=2)
# Show the predicted expectation from BiGQD.density():
points(res$cumulants[5,i]~res$cumulants[1,i],bg='white',pch=21,cex=1.5)
axis(1);axis(2);
# Add a legend:
legend('topright',lty=c('dashed',NA),pch=c(NA,21),lwd=c(2,NA),
legend=c('Simulated Expectation','Predicted Expectation'))
})
}


## ----eval=FALSE----------------------------------------------------------
#  browseVignettes('DiffusionRgqd')

## ----eval=FALSE----------------------------------------------------------
#  demo(package = 'DiffusionRgqd')

