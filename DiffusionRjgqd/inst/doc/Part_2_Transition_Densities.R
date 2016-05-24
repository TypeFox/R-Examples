## ------------------------------------------------------------------------
library(DiffusionRjgqd)
# Some parameter values:
a <- 1; b <- 5; sigma <- 0.25; lam_0 <- 0.5; nu <- 0.3; 

# Define the model:

# Diffuse part
G0 <- function(t){a*b}
G1 <- function(t){-a}
Q1 <- function(t){sigma}
# Jump part
Lam0  <- function(t){lam_0}
Jlam  <- function(t){nu}

## ----fig.align = 'center'------------------------------------------------
 Xt <- seq(2,9,1/20)
 Xs <- 6
 t  <- 0.4
 s  <- 0
 dt <- 1/100

 M <- JGQD.density(Xs,Xt,s,t,dt, Jdist = 'Exponential', factorize = TRUE)

 persp(x=M$Xt,y=M$time,z=M$density,col='green',xlab='State (X_t)',ylab='Time (t)',
      zlab='Density f(X_t|X_s)',border=NA,shade=0.5,theta=145)

## ------------------------------------------------------------------------
library(DiffusionRjgqd)
# Some parameter values:
kap <- 1; beta_0 <- 5; beta_1 <- 3; sigma <- 0.15; 
lam_0 <- 0.5; lam_1 <- 0.1; mu_z <- 0.5; sigma_z <- 0.2; 

# Define the model:
JGQD.remove()
# Diffuse part
G0 <- function(t){kap*(beta_0+beta_1*sin(2*pi*t))}
G1 <- function(t){-kap}
Q1 <- function(t){sigma^2}
# Jump part
Lam0  <- function(t){lam_0}
Lam1  <- function(t){lam_1}
Jmu   <- function(t){mu_z}
Jsig  <- function(t){sigma_z}

## ----fig.align = 'center'------------------------------------------------
 Xt <- seq(2,9,1/20)
 Xs <- 6
 t  <- 5
 s  <- 0
 dt <- 1/100

 M <- JGQD.density(Xs,Xt,s,t,dt, Jdist = 'Normal')

 persp(x=M$Xt,y=M$time,z=M$density,col='green',xlab='State (X_t)',ylab='Time (t)', zlab='Density f(X_t|X_s)',border=NA,shade=0.5,theta=145)

## ----fig.align = 'center'------------------------------------------------
 Xt <- seq(5.5,8,1/50)
 t  <- 0.2
 M2  <- JGQD.density(Xs,Xt,s,t,dt, Jdist = 'Normal',factorize = TRUE, print.output = FALSE)
 
 # Plot the transitional density at t =0.2 and the evolution of the zero
 # jump probability
 par(mfrow=c(1,2))
 plot(M2$density[,20]~M2$Xt,type='l',main ='Transition density at t = 0.2')
 plot(M$zero_jump_prob~M$time,type='l',main =expression(P(N[t]-N[s]==0)))
 
 # Superimpose the short horizon on the probability trajectory:
 abline(v=t,lty='dotted')

## ----fig.align= 'center'-------------------------------------------------
 # Package for purely diffuse GQDs:
 library(DiffusionRgqd)
 
 M3  <- GQD.density(Xs,Xt,s,t,dt, print.output = FALSE)
 
 plot(M2$density[,20]~M2$Xt,type='l',ylim=c(0,2.8),main ='Transition density at t = 0.2')
 lines(M3$density[,20]~M3$Xt,col='blue',lty='dashed')
 legend('topright',legend = c('Jump','Diffuse'),col=c(1,4),lty=c('solid','dashed'))
 

## ----fig.align='center'--------------------------------------------------
JGQD.remove()

G0=function(t){2*5+2*sin(1*pi*t)}
G1=function(t){-2}
Q1=function(t){1}

l1 <- 1
l2 <- 3
rho1 <- 0.25
rho2 <- 1

Jmu  <- function(t){1}
Jsig <- function(t){0.25}
Lam0 <- function(t){l1*(rho2+rho1*exp(-(rho1+rho2)*t))/(rho1+rho2)+l2*rho1/(rho1+rho2)*(1-exp(-(rho1+rho2)*t))}
t <- seq(0,5,1/100)
# Intensity assumes the expectation trajectory of the intensity process:
plot(Lam0(t)~t,type='l')

TT  <- 5
res <- JGQD.density(4,seq(2,14,1/10),0,TT,1/100)

## ------------------------------------------------------------------------
#' Now simulate the jump diffusion

mu     <- function(x,t){G0(t)+G1(t)*x} # Drift
sigma  <- function(x,t){sqrt(Q1(t)*x)} # Diffusion
j      <- function(x,z){z}             # Jumps
simulate <- function(x0=4,N=25000,pts = 1:4)
{
  d=0
  delta=1/1000
  tt=seq(0,TT,delta)
  X=rep(x0,N)

  kkk=1
  MM= matrix(0,4,length(tt))
  MM[,1]=X[1]^{1:4}

  xtrak = rep(x0,length(tt))
  jtrak = rep(0,length(tt))
  etrak = rep(0,length(tt))
  ltrak = rep(0,length(tt))
  sttes = rep(l1,N)
  L=list()
  kkk =1
  p.states = jtrak
  e.states = jtrak
  for(i in 2:length(tt))
  {
    X=X+mu(X,d)*delta+sigma(X,d)*rnorm(N,sd=sqrt(delta))
    d=d+delta

    events = ((sttes)*delta*exp(-(sttes)*delta)>runif(N))
    xpre=X[1]
    if(any(events))
    {
      wh=which(events)
      X[wh]=X[wh]+j(X[wh],rnorm(length(wh),Jmu(d),Jsig(d)))

    }
    jtrak[i] = X[1]-xpre
    xtrak[i]  = X[1]


    prbs1=rho1/(rho1+rho2)*(1-exp(-(rho1+rho2)*delta))
    prbs2=rho2/(rho1+rho2)*(1-exp(-(rho1+rho2)*delta))

    whh1=which(sttes==l1)
    whh2=which(sttes==l2)
    if(length(whh1)!=0)
    {
       whh1.2=which(runif(length(whh1))<prbs1)
       if(length(whh1.2)!=0)
       {
        sttes[whh1][whh1.2]=l2
       }
    }
    if(length(whh2)!=0)
    {
      whh2.2=which(runif(length(whh2))<prbs2)
      if(length(whh2.2)!=0)
      {
       sttes[whh2][whh2.2]=l1
      }
    }

    p.states[i] = mean(sttes==l2)
    e.states[i] = mean(sttes)
    etrak[i]   = sttes[1]

    MM[1,i]=sum(X)/N
    MM[2,i]=sum(X^2)/N
    MM[3,i]=sum(X^3)/N
    MM[4,i]=sum(X^4)/N
    if(sum(round(pts,3)==round(d,3))!=0)
    {
       L[[kkk]] = hist(X,plot=F,breaks=25)
       kkk=kkk+1
    }
  }
  return(list(MM=MM,tt=tt,X=X,xtrak=xtrak,etrak=etrak,jtrak=jtrak,ltrak=ltrak,hists=L,p.states=p.states,e.states=e.states))
}
 res2 <- simulate()

## ----fig.align='center'--------------------------------------------------
 par(mfrow=c(3,1))
 plot(res2$xtrak~res2$tt,type='l',col='blue',main='Trajectory',xlab = 'time',ylab ='X_t')
 plot(res2$jtrak~res2$tt,type='h',col='black',ylim=c(0,3),lwd=2,main='Jumps',xlab ='time',ylab ='Z_t')
 plot(res2$etrak~res2$tt,type='s',col='black',ylim=c(0,6),lwd=1,main='Intensities',xlab ='time',ylab ='Z_t')
 abline(h =c(l1,l2),lty='dotted',col='lightgrey')

 
 par(mfrow=c(2,2))
 for(i in 1:4)
 {
  plot(res2$MM[i,]~res2$tt,type='l',main='Moment trajectory',xlab='Time (t)',
       ylab=paste0('m_',i,'(t)'))
  lines(res$moments[i,]~res$time,lty='dashed',col='red',lwd=2)
 }

## ----fig.align='center'--------------------------------------------------
 persp(x=res$Xt,y=res$time,z=res$density,col='green',xlab='State (X_t)',ylab='Time (t)', zlab='Density f(X_t|X_s)',border=NA,shade=0.5,theta=145)
 
 par(mfrow=c(2,2))
 for(i in 1:4)
 {
   plot(res2$hists[[i]]$density~c(res2$hists[[i]]$mids-diff(res2$hists[[i]]$mids)[1] / 2),
        type = 's',lty = 'solid', lwd = 1, xlab = 'time', 
        ylab = 'Density', main = paste('Density at time t =',i))
   lines(res$density[,i*100]~res$Xt,col='darkblue')
 }

## ----fig.align='center'--------------------------------------------------
# Shrink the diffusion coefficient:
Q1=function(t){0.2}
 
# Approximate, but use factorization:
TT  <- 0.5
res <- JGQD.density(4,seq(2,8,1/10),0,TT,1/100, factorize = TRUE)

# Re-simulate and record histograms at new points:
res2 <- simulate(pts=seq(0.1,0.4,0.1))

persp(x=res$Xt,y=res$time,z=res$density,col='green',xlab='State (X_t)',ylab='Time (t)', zlab='Density f(X_t|X_s)',border=NA,shade=0.5,theta=145)

par(mfrow=c(2,2))
for(i in 1:4)
{
  plot(res2$hists[[i]]$density~c(res2$hists[[i]]$mids-diff(res2$hists[[i]]$mids)[1] / 2),
       type = 's',lty = 'solid', lwd = 1, xlab = 'time', 
       ylab = 'Density', main = paste('Density at time t =',i*0.05))
  lines(res$density[,i*10]~res$Xt,col='darkblue')
}

## ----eval=FALSE----------------------------------------------------------
#  browseVignettes('DiffusionRjgqd')

