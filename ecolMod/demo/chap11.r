##############################
## Soetaert and Herman      ##
## ecological modelling     ##
## Figures from chapter 11  ##
## Testing and Validating   ##
##############################

opar <- par()
par(ask=TRUE)
par(mfrow=c(1,1))
figNr <- 1
subtitle <- function()
{
 mtext(side=1,outer=TRUE,"Soetaert and Herman - chapter 11  ",cex=0.7,adj=1,line=-1.5)
 mtext(side=1,outer=TRUE,paste("  Fig. 11.",figNr,sep=""),cex=0.7,adj=0,line=-1.5)
 figNr <<- figNr +1
}

###############################################################################
####======================================================================#####
####                                Theory                                #####
####======================================================================#####
###############################################################################

########################################
# Figure 11.1. BOD model scheme
########################################

par(mar=c(0,0,0,0))
par(mfrow=c(1,1))
emptyplot()
rect(0.1,0.1,0.9,0.9)
rect(0.11,0.11,0.89,0.6,border=NA,col="lightgrey")
elpos <- matrix(nrow=2,data=c(0.25,0.75,0.4,0.4))
straightarrow(from=elpos[1,]+c(0,0.1),to=elpos[1,]+c(0,0.3),arr.pos=1)
straightarrow(to=elpos[1,]+c(0,0.1),from=elpos[1,]+c(0,0.3),arr.pos=1,arr.adj=1)
text(elpos[1,1],elpos[1,2]+0.4,expression(O[2]^"*"))
bentarrow (from=elpos[1,],to=elpos[1,]+c(0.2,-0.2),arr.pos=1.0,lwd=2)
bentarrow (from=elpos[2,],to=elpos[2,]+c(-0.2,-0.2),arr.pos=1.0,lwd=2)

textrect(elpos[1,],0.1,0.1,lab=expression(O[2]),shadow.size=0.02)
textrect(elpos[2,],0.1,0.1,lab="BOD",shadow.size=0.02)
subtitle()

########################################
# Figure 11.2. BOD model
########################################

par(oma=c(0,0,2,0))
par(mar=c(5.1,4.1,4.1,2.1))
par(mfrow=c(2,2))
k       = 0.1             # /day     - reaeration
O2sat   = 300             # mmol/m3  - saturated oxygen concentration
r       = 0.05            # /day     - BOD decay rate
O2_0    = 250             # mmol/m3  - Initial oxygen concentration
BOD_0   = 500             # mmol/m3  - Initial BOD concentration
ks      = 0               # mmol/m3  - half-saturation concentration

# numerical model
numBOD <- function (time,state,pars)
{
 with (as.list(state), 
  {
    dO2  <- -r*BOD*O2/(O2+ks)+k*(O2sat-O2)
    dBOD <- -r*BOD*O2/(O2+ks)
    return(list(c(dO2,dBOD)))
  }
      )
}

# analytical solution for O2
analytical <- function(x,k=0.1,r=0.05,O2sat=300) 
 BOD_0*r*(exp(-k*x)-exp(-r*x)) /(k-r)+O2_0*exp(-k*x)+O2sat*(1-exp(-k*x))
 
# A comparison numerical / analytical model
# numerical solution plotted as points
times <- 0:100
state <- c(O2=O2_0,BOD=BOD_0)
out   <- as.data.frame(ode(state,times,numBOD,0))
plot(out$time,out$O2,xlab="time",ylab="mmol O2/m3",lwd=2,main="Correctness of solution") 

# analytical solution - added as a curve
curve(analytical(x,k),lty=1,lwd=1,add=TRUE)
legend("bottomright",c("analytic","numerical"),lwd=c(1,2),lty=c(1,NA),pch=c(NA,1))
writelabel("A")

# B: internal logic
# wrong use of model : too low reaeration -> negative concentration

k     <- 0.01
times <- 0:200
state <- c(O2=O2_0,BOD=BOD_0)
out   <- as.data.frame(ode(state,times,numBOD,0))
plot(out$time,out$O2,xlab="time",ylab="mmol O2/m3",main="Internal logic",type="l",lty=2) 
abline(h=0,lty=3)

ks    <- 1
state <- c(O2=O2_0,BOD=BOD_0)
out2  <- as.data.frame(ode(state,times,numBOD,0))
lines(out2$time,out2$O2,lwd=2)
legend("bottomright",c("no O2 limitation","O2 limitation"),lwd=c(1,2),lty=c(2,1))
writelabel("B")

# C: global sensitivity
k       <- 0.1           
rseq   <- seq(0.0,0.2,by=0.002)
rseq   <- rseq[rseq!=k]    # cannot calculate analytical solution for this...
minO2  <- rseq
for (i in 1:length(rseq)) 
  minO2[i]  <- min(analytical(times,r=rseq[i]))
plot(rseq,minO2,type="l",lwd=2 ,xlab="r, /day",ylab="minimum O2, mmol/m3",main="global sensitivity")
writelabel("C")

mtext(side=3,outer=TRUE,line=0,"BOD-O2 model",cex=1.25,font=2)

# D: local sensitivity

times   <- 0:100 
ss      <- 1.1

kp      <- k * ss         # /day     - reaeration
O2satp  <- O2sat*ss        # mmol/m3  - saturated oxygen concentration
rp      <- r*ss           # /day     - BOD decay rate

ref  <- analytical(times)
outk <- analytical(times,k=kp)
outs <- analytical(times,O2sat=O2satp)
outr <- analytical(times,r=rp)
outm <- mean(ref)
ss   <- cbind(k=(outk-ref)/outm/0.1,sat=(outs-ref)/outm/0.1,r=(outr-ref)/outm/0.1)

plot(times,ref,ylim=range(c(ref,outs)),type="l",lwd=2,xlab="time",ylab="mmol O2/m3",main="local sensitivity")
lines(times,outs,lwd=2,lty=2)
arrseq <- seq(10,100,10)#c(10,30,50,70,90)
Arrows(times[arrseq],ref[arrseq],times[arrseq],outs[arrseq],arr.len=0.25,arr.adj=1)
legend("topleft",c(expression(O[2]^"*"== 300),expression(O[2]^"*"== 330)),lwd=2,lty=c(1,2))

writelabel("D")

par(new=TRUE)
par(fig=c(0.7,0.99,0.01,0.35))                                
plot(times,ss[,2],type="l",lwd=2,  
   xlab="",ylab="",axes=FALSE,frame.plot=TRUE)
points(times[arrseq],ss[arrseq,2])
text(mean(times),diff(range(ss[,2]))/2,expression(S["i,j"]))
#

msqr <- sqrt(colSums(ss*ss)/length(times))

par(fig=c(0,1,0,1))
subtitle()


###############################################################################
####======================================================================#####
####                            R case studies                            #####
####======================================================================#####
###############################################################################


########################################
# Fig 11. 3. conceptual model O2 budget
########################################

par(mar=c(2,2,2,2))
emptyplot(c(-0.2,0.2),c(-0.5,0.5),asp=FALSE,main="")
mid <- c(-0.1,0)

rseq1 <- seq(0.3,1,length.out=50)
col  <- grey(c(rseq1,rev(rseq1)))
rx<-0.05
ry<-0.1
filledcylinder(rx=rx,ry=ry,len=0.85,angle=90,col=col,mid=mid,lcol="black",lwd=1)
for (i in seq(0.2,0.8,by=0.2))plotellipse(mid=c(mid[1],0.425),rx=i*ry,ry=i*rx,lwd=1,lcol="grey")
straightarrow (from=c(mid[1]+0.17,0.0),to=mid,lwd=2,arr.pos=0.5,arr.length=0.5)
bentarrow (to=c(mid[1]-0.8*ry,-0.2),from=mid,lwd=2,arr.pos=1,arr.length=0.3,arr.type="triangle")
textrect   (mid,0.05,0.05,lab=expression(O[2]),shadow.size=0.005,cex=1.5)

text( mid[1]-0.6*rx,-0.2,"Q",cex=1.15)
text( mid[1]+0.21,0,"J",cex=1.15)
text(0.125,0.35,expression(BWO[2]))
par(new=TRUE)
par(fig=c(0.45,1,0.5,0.85))
plot(x<-seq(0,2*pi,len=20),sin(x),xlab="",ylab="",axes=FALSE,type="l",
     xaxs="i",yaxs="i",lwd=2,cex.main=1)
box(col="grey")
subtitle()

par(fig=c(0,1,0,1))
BW     <- 2            # mmol/m3,       oxygen concentration in surrounding water
Da     <- 0.5          # cm2/d          effective diffusion coefficient in organism
R      <- 0.0025        # cm             radius of organism
Q      <- 250000       # nM/cm3/d       oxygen consumption rate per volume per day
L      <- 0.05          # cm             length of organism

# the analytical cylindrical model
cylinder <- function(Da,Q,BW,R,r)  BW+Q/(4*Da)*(r^2-R^2)

par(mar=c(5.1,4.1,4.1,2.1))
par(mfrow=c(2,2))

N  <- 10                           # we consider 10 layers in the nematode body
dx <- R/N                          # thickness of each layer
x  <- seq(dx/2,by=dx,length.out=N) # distance of center to mid-layer
xi <- seq(0,by=dx,length.out=N+1)  # distance to layer interface
dxi<- c(rep(dx,N),dx/2)            # dispersion distances
A   <- 2*pi*x *L                   # surface at mid-layer depth
Ai  <- 2*pi*xi*L                   # surface at layer interface

# the numerical model

oxygen <- function (time, O2, pars)
#==============================================================================
# the rate of change of oxygen in the body of a cylindrical organism
# outer boundary = concentration boundary
# lower boundary = zero-gradient boundary
#==============================================================================

  {
    BWO2 <- BW*(1-0.8*sin(2*pi*time*24))
    # outer concentration imposed (BW), lower: zero gradient
    Flux    <- -Da * diff(c(O2[1],O2,BWO2))/dxi   #diffusion

    # Rate of change = Flux gradient - oxygen consumption
    dO2     <- -diff(Ai * Flux)/A/dx -Q

    return (list(dO2=dO2,c(Flux=Flux,BWO2=BWO2)))
  }


CONC  <- steady.1D (runif(N),func=oxygen,nspec=1,atol=1e-10,time=0)
O2    <- CONC$y
plot(x,O2,xlab="distance, cm",ylab="oxygen, mumol/l")
lines(x, BW+Q/(4*Da)*(x^2-R^2))
legend ("topleft",lty=c(1,NA),pch=c(NA,1),c("analytical solution","numerical approximation"))
writelabel("A")

N  <- 100                          # we consider 100 layers in the nematode body
dx <- R/N                          # thickness of each layer
x  <- seq(dx/2,by=dx,length.out=N) # distance of center to mid-layer
xi <- seq(0,by=dx,length.out=N+1)  # distance to layer interface
dxi<- c(rep(dx,N),dx/2)            # dispersion distances
A   <- 2*pi*x *L                   # surface at mid-layer depth
Ai  <- 2*pi*xi*L                   # surface at layer interface

CONC  <- steady.1D (runif(N),func=oxygen,nspec=1,atol=1e-10,time=0)
O2    <- CONC$y
plot(x,O2,xlab="distance, cm",ylab="oxygen, mumol/l")
lines(x, BW+Q/(4*Da)*(x^2-R^2))
legend ("topleft",lty=c(1,NA),pch=c(NA,1),c("analytical solution","numerical approximation"))
writelabel("B")


times<- seq(0,1/24,length.out=120)
out   <-as.data.frame(ode.band(O2,times,func=oxygen,parms=0,rtol=1e-10,atol=1e-10,nspec=1))
oxy   <- out[,2:101]
plot(times*24,out$BWO2,xlab="time, hour",ylab="mumol/l", main="BW concentration",type="l", lwd=2)
writelabel("C")
image(times*24,y=x,z=as.matrix(oxy),xlab="time, hour",ylab="distance ",main="Dynamic simulation",col=femmecol(100))
contour(times*24,y=x,z=as.matrix(oxy),add=TRUE)
box()
writelabel("D")
mtext(outer=TRUE,side=3,"Oxygen in cylindrical body",cex=1.5,line=-1.5)
subtitle()

########################################
# Fig. 11.5 BACT GLUCOSE
########################################

par(mfrow=c(2,2))
par(mar=c(2,2,2,2))
emptyplot()
rect(0.05,0.05,0.95,0.95)
rect(0.06,0.06,0.94,0.75,border=NA,col="lightgrey")
elpos <- matrix(nrow=2,data=c(0.3,0.7,0.4,0.4))
dar <- c(0,0.025)
straightarrow (from=elpos[1,]+dar,to=elpos[2,]+dar,arr.pos=0.6,lwd=2)
straightarrow (to=elpos[1,]-dar,from=elpos[2,]-dar,arr.pos=0.6,lwd=2)
bentarrow(elpos[1,],elpos[1,]+c(-0.15,-0.15),arr.type="circle",arr.pos=1,arr.length=0.3)
textrect(elpos[1,],0.1,0.1,lab="BACT",shadow.size=0.02)
textrect(elpos[2,],0.1,0.1,lab="GLUC",shadow.size=0.02)
par(mar= c(5.1,4.1,4.1,2.1))
writelabel("A",at=0.15,line=0.9)


pars <- list(Bini=0.1,Sini=100,gmax =0.5,eff = 0.5,
              ks =0.5, rB =0.01, dB =0.01)
model <- function(t,state,pars)
{
with (as.list(c(state,pars)), {
dBact = gmax*eff*Sub/(Sub+ks)*Bact - dB*Bact - rB*Bact
dSub  =-gmax    *Sub/(Sub+ks)*Bact + dB*Bact
return(list(c(dBact,dSub)))
                             })
}

tout    <- seq(0,50,by=0.5)
state   <- c(Bact=pars$Bini,Sub =pars$Sini)
out     <- as.data.frame(ode(state,tout,model,pars))
plot(out$time,out$Bact,ylim=range(c(out$Bact,out$Sub)),
     xlab="time, hour",ylab="molC/m3",type="l",lwd=2)
lines(out$time,out$Sub,lty=2,lwd=2)
lines(out$time,out$Sub+out$Bact)
legend("topright",c("Bacteria","Glucose","TOC"),
       lty=c(1,2,1),lwd=c(2,2,1))
writelabel("B",at=-0.1,line=0.9)

subtitle()

par(mfrow=c(1,1))


###############################
# Fig. 11.6 graphical analysis 
###############################

Reference <- out
pp        <- unlist(pars)
tiny      <- 0.1
par(mfrow=c(3,3))
for (i in 1:length(pars))
{
 pars[i] <- pp[i]*(1+tiny)
 state   <- c(Bact=pars$Bini,Sub =pars$Sini)
 out     <- as.data.frame(ode(state,tout,model,pars))

 plot(out$time,out$Bact,xlab="hour",ylab="molC/m3",
      type="l",lwd=2,main=names(pars)[i])
 lines(Reference$time,Reference$Bact,lty=2)
 pars[i] <- pp[i]
}
plot(0,axes=FALSE,xlab="",ylab="",type="n")
legend("center",c("perturbed","reference"),lwd=c(2,1),lty=c(1,2))
subtitle()

###############################
# Fig. 11.7 Sensitivity functions
###############################

state <- c(Bact=pars$Bini, Sub =pars$Sini)
yRef  <- as.data.frame(ode(state,tout,model,pars))$Bact
pp    <- unlist(pars)
nout  <- length(yRef)
npar    <- length(pars)

# perturbation factor
tiny    <- 1e-8
dp      <- pp*tiny

Sens    <- matrix(nrow=nout,ncol=npar,NA)

for (i in 1:npar)
{
  dval    <- pp[i]+dp[i]
  pars[i] <- dval
  state   <- c(Bact=pars$Bini,Sub =pars$Sini)
  yPert   <- as.data.frame(ode(state,tout,model,pars))$Bact
  Sens[,i]<- (yPert-yRef)/tiny
  pars[i] <- pp[i]
}
colnames(Sens) <- names(pars)
rownames(Sens) <- tout
format(as.data.frame(Sens[1:5,]),digits=2)

par(mfrow=c(1,1))
matplot(tout,Sens,type="l",lty=1 :10,col=1:10)
legend("topright",names(pars),lty=1:10,col=1:10)
subtitle()

mabs <- colMeans(abs(Sens))
msqr <- sqrt(colSums(Sens*Sens)/nout)
format(data.frame(msqr,mabs),digits=2)

###############################
# Fig. 11.8 bivariate analysis 
###############################

par(mar=c(5.1,4.1,4.1,2.1))
par(oma=c(0,0,0,0))
panel.cor <- function(x, y)
             text(x=mean(range(x)),y=mean(range(y)),
             labels=format(cor(x,y),digits=2))

pairs(Sens,upper.panel=panel.cor)
mtext(outer=TRUE,side=3,line=-1,
      "Sensitivity functions",cex=1.5)
subtitle()

par(ask=opar$ask)
par(mar=opar$mar)
par(oma=opar$oma)
par(mfrow=opar$mfrow)

