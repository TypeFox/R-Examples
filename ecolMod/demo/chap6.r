##############################
## Soetaert and Herman      ##
## ecological modelling     ##
## Figures from chapter 6   ##
## Numerical solutions      ##
##############################

opar <- par()
par(ask=TRUE)
par(mfrow=c(1,1))
figNr <- 1
subtitle <- function()
{
 mtext(side=1,outer=TRUE,"Soetaert and Herman - chapter 6  ",cex=0.7,adj=1,line=-1.5)
 mtext(side=1,outer=TRUE,paste("  Fig. 6.",figNr,sep=""),cex=0.7,adj=0,line=-1.5)
 figNr <<- figNr +1
}

###############################################################################
####======================================================================#####
####                                Theory                                #####
####======================================================================#####
###############################################################################

#####################################################
# Fig. 6.1. numerical approximation in time and in space 
#####################################################
par(mfrow=c(1,1))
par(mar=c(0,0,0,0))
emptyplot(xlim = c(-0.25,1.25),ylim = c(-0.1,0.7))

basemid <-  c(0,0.1)
dy <- 0.1
ds <- 0.05    ;r2=0.15
col1 <- rev(grey(seq(0.4,0.6,len=4)))
col2 <- rev(grey(seq(0.6,0.8,len=4)))
col3 <- rev(grey(seq(0.8,1.0,len=4)))

for (i in 1:4) {
filledcylinder(rx=ds,ry=r2,len=dy, mid=basemid+c(0.,(i-1)*dy),angle=90,dr=0.01,
             col=col1[i],lwd=2,topcol=col1[i],lcol="black")
                }

basemid <-  c(0.5,0.1)
for (i in 1:4) {
filledcylinder(rx=ds,ry=r2,len=dy, mid=basemid+c(0.,(i-1)*dy),angle=90,dr=0.01,
             col=col2[i],lwd=2,topcol=col2[i],lcol="black")
                }

basemid <-  c(1.,0.1)

for (i in 1:4) {
filledcylinder(rx=ds,ry=r2,len=dy, mid=basemid+c(0.,(i-1)*dy),angle=90,dr=0.01,
             col=col3[i],lwd=2,topcol=col3[i],lcol="black")
                }

text(0,0.6,expression(t),cex=1.5,font=2)
text(0.5,0.6,expression(t+Delta~t),cex=1.5,font=2)
text(1.,0.6,expression(t+2*Delta~t),cex=1.5,font=2)

Arrows(0.1,0.6,0.2,0.6,arr.type="triangle",arr.length=0.15,segment=TRUE)
Arrows(0.68,0.6,0.78,0.6,arr.type="triangle",arr.length=0.15,segment=TRUE)
subtitle()
#####################################################
# Fig. 6.2.  Taylor expansion                          
#####################################################


par(mfrow= c(2,2))
par(mar=c(5.1,4.1,4.1,2.1))
# function and derivates
fun <- function (x) 2*x^3
der <- function (x) 6*x^2  # 1st derivate
der2<- function (x) 12*x   # 2nd derivate
der3<- function (x) 12     # 3rd derivate

# points
x1  <- 1
hh  <- 1
x2  <- x1+hh
y1  <- fun(x1)
y2  <- fun(x2)

# main plot function
plotf <- function()
{
curve(fun(x),0.5,2.5,axes=FALSE,xlab="",ylab="",lwd=2,ylim=c(-2,32))
axis(1,labels=FALSE,at=c(0,3))
axis(1,at=c(1,2),labels=c("x","x+h"),font=2,cex.axis=0.8)
axis(2,labels=FALSE,at=c(-5,35))
axis(2,at=c(y1,y2),labels=c("f(x)","f(x+h)"),font=2,cex.axis=0.8)
points(x1,y1,pch=16,col="black",cex=2)
points(x2,y2,pch=16,col="black",cex=2)
segments(x1,-5,x1,y1,lty=2)
segments(x1,y1,-5,y1,lty=2)
segments(x2,-5,x2,ye,lty=2,col="grey")
segments(x2,y2,-5,y2,lty=2)
if (ye != y1) segments(x2,ye,-5,ye,lty=2,col="grey")
}

# 0-th order expansion
ye  <- fun(x1)
plotf()
segments(x1,y1,x2,ye,lty=1)
points(x2,ye,pch=16,col="darkgrey",cex=2)
text(1.5,30,expression(f[x+h]%~~%f[x]),cex=1.25)
writelabel("A")

# 1-st order expansion
ye <- fun(x1) + hh*der(x1)
plotf()
hhs <- seq(0,hh,len=100) 
xx  <- x1+hhs
lines(xx,fun(x1)+hhs*der(x1))
points(x2,ye,pch=16,col="darkgrey",cex=2)
text(1.5,30,expression(f[x+h]%~~%f[x]+h*f[x]^" '") ,cex=1.15)
writelabel("B")

# 2-nd order expansion
ye <- fun(x1) + hh*der(x1) + hh^2/2* der2(x1)
plotf()
hhs <- seq(0,hh,len=100)

lines(xx,fun(x1)+hhs*der(x1)+hhs^2/2*der2(x1))
points(x2,ye,pch=16,col="darkgrey",cex=2)
text(1.5,28,expression(f[x+h]%~~%f[x]+h*f[x]^" '"+frac(h^2,2)*f[x]^" ''") ,cex=1.15)
writelabel("C")

# 3-rd order expansion
ye <- fun(x1) + hh*der(x1) + hh^2/2* der2(x1)   + hh^3/6* der3(x1)
plotf()
hhs <- seq(0,hh,len=100)
lines(xx,fun(x1)+hhs*der(x1)+hhs^2/2*der2(x1)+ hhs^3/6* der3(x1))
points(x2,ye,pch=16,col="darkgrey",cex=2)
text(1.5,28,expression(f[x+h]%~~%f[x]+h*f[x]^" '"+frac(h^2,2)*f[x]^" ''"+frac(h^3,3~"!")*f[x]^" '''") ,cex=1.15)
writelabel("D")
subtitle()
#####################################################
# Fig. 6.3. numerical errors 
#####################################################

par(mfrow=c(2,2),cex.axis=0.8,cex.lab=0.9)

mu <-0.5
curve(1*exp(mu*x),0,10,xlab="time",ylab="y" )
x0<- 1
dt<-0.1
times<- seq(0,10,by=dt)
x <- times
x[1]<-x0
for (i in 2:length(times))x[i]<-x[i-1]+dt*mu*x[i-1]
points(times,x)
legend("topleft",legend=c("true solution","approximation"),
cex=0.9,lty=c(1,NA),pch=c(NA,1))

writelabel("A")


pMax<-1
ks  <-1
Nin<-10
r  <- 0.24
pini <- 1
nini <- 0.1

model<-function(p,n)
{  #estimates the rate of change
dp <- pMax*n/(n+ks)*p-r*p
dn <- -pMax*n/(n+ks)*p+r*(Nin-n)
return(list(dp=dp,dn=dn))
}

euler <- function(dt)
{
times<- seq(0,20,by=dt)
PP <- NN <- times
PP[1]<-pini
NN[1]<-nini
for (i in 2:length(times))
{
rate <-model(PP[i-1],NN[i-1])
PP[i]<-PP[i-1]+dt*rate$dp
NN[i]<-NN[i-1]+dt*rate$dn
}
return(list(TT=times,PP=PP,NN=NN))
}

EU1<-euler(0.01)
EU2<-euler(0.5)

plot(EU1$TT,EU1$NN,type="l",ylim=range(c(EU1$NN,EU2$NN)),xlab="time" ,ylab="y" )
abline(h=0,lty=2)
points(EU2$TT,EU2$NN,type="b" ,pch=21,bg=c("white","black")[(EU2$NN<0)+1],col="black")
legend("topright",legend=c("true solution","approximation"),cex=0.9,lty=1,pch=c(NA,1))
writelabel("B")
subtitle()

#####################################################
# Fig. 6.4. numerical integration scheme
#####################################################

par(mfrow=c(1,1))
par(mar=c(0,0,0,0))

openplotmat()
elpos <-coordinates (c(1,1,1,1,1))
segmentarrow(elpos[4,],elpos[2,],arr.side=3,arr.pos=0.15,dd=0.4)

dy  <- 0.0175
dx  <- 0.15

for (i in 1:4) straightarrow (elpos[i,],elpos[i+1,],lwd=2,arr.pos=0.5)
textround  (elpos[1,],dx,3*dy,lab=c("Initialise model","t=0",expression(C=C[0])))
textround  (elpos[2,],dx,3*dy,lab=c("Rate of change","",expression(over(dC,dt)==sources-sinks),""))
textround  (elpos[3,],dx,3*dy,lab=c("Integrate",expression(C[t]%=>%C[t+Delta~t]),expression(t %=>%t+Delta~t)))
textdiamond(elpos[4,],0.085,lab="t=tend?")
textround  (elpos[5,],dx,dy*1.2 ,lab="Output")

text(0.2 ,elpos[4,2]+0.015,"No",adj=0,font=3)
text(0.52,elpos[4,2]-0.07 ,"yes",adj=0,font=3)
subtitle()

#####################################################
# Fig. 6.5. numerical error of euler 
# and runge-kutta integration method   
#####################################################


par(mfcol=c(2,2))
par(mar=c(3,6,3,0),las=1)

mu <- 0.9
fun <- function(x) exp(mu*x)
der <- function(x) mu*exp(mu*x)

curve(fun(x),0,2.8,lwd=2,axes=FALSE,frame.plot=TRUE,xlab="",ylab="",
       ylim=c(0,8),main="Explicit euler method")

dx  <- 1
x1  <- 1
x2  <- x1+dx
y1  <- fun(x1)
dx1 <- der(x1)

segments(x1,0,x1,y1,lty=2)
y2  <- y1 +dx*dx1
segments(x1,y1,x2,y2,lty=1)
segments(x2,0,x2,y2,lty=2)

y2t  <- fun(x2)
segments(x1,y1,x2,y2,lty=1)
segments(x2,0,x2,y2t,lty=2)
segments(x2,y1,x1,y1,lty=2)

points(x1,y1,pch=16,col="black",cex=2)
points(x2,y2,pch=16,col="darkgrey",cex=2)
points(x2,y2t,pch=18,col="black",cex=2)

axis(2,at=c(y1,y2,y2t),labels=c("f(t1)",expression(f(t1)+Delta~t*over(df,dt)(t1)),"f(t2)"),font=2,cex.axis=0.8)
axis(1,at=c(x1,x2),labels=c("t1","t2"))
text(x2+0.05,y2,adj=0,"estimated",col="darkgrey")
text(x2+0.05,y2t,adj=0,"true value")
text(0.6*x1+0.4*x2,mean(c(y1,y2))-0.5,adj=0,srt=40,"tangent line",cex=0.8)
text(mean(c(x1,x2)),y1+0.25,adj=0,expression(Delta~t))
writelabel("A")


# implicit Euler method
mu <- 0.9
fun <- function(x) exp(mu*x)
der <- function(x) mu*exp(mu*x)

dx  <- 1
x1  <- 1
x2  <- x1+dx

y1  <- fun(x1)
y2  <- fun(x2)

yse  <- y1 + dx*0.5*(der(x1) +der(x2))          # semi-implicit
yim <- y1 + dx* der(x2)                         # implicit

y1  <- fun(x1)
dx1 <- der(x1)
dx2 <- der(x2)

curve(fun(x),0,2.8,lwd=2,axes=FALSE,frame.plot=TRUE,xlab="",ylab="",
       ylim=c(0,8),main="Implicit euler method")
segments(x1,0,x1,y1,lty=2)
axis(2,at=c(y1,yim,y2),labels=c("f(t1)",expression(f(t1)+Delta~t*over(df,dt)(t2)),"f(t2)"),font=2,cex.axis=0.8)
axis(1,at=c(x1,x2),labels=c("t1","t2"))

segments(x1,y1,x2,yim,lty=1)
segments(x2,0,x2,yim,lty=2)
dd <- 0.5
segments(x2-dd,y2-dd*dx2,x2+dd,y2+dd*dx2,lty=1)
#segments(x1-dd,y1-dd*dx1,x1+dd,y1+dd*dx1,lty=2)

points(x1,y1,pch=16,col="black",cex=2)
points(x2,yim,pch=16,col="darkgrey",cex=2)
points(x2,y2,pch=18,col="black",cex=2)
writelabel("C")

#=================================================================
## SIMPLE RUNGEKUTTA INTEGRATOR - NO ERROR CHECKING
#=================================================================

NewY <- function(y,t,derivative,dt)

#----------------------------------------------------------------
# One runge-kutta step RK4, performed at time t
# y:          value of state variable at time t
# derivative: returns the rate of change of y, as a VECTOR 
# dt:         the fixed time step
# Returns the updated value of the state variables
#----------------------------------------------------------------
     {
        # four steps
        y1 <- dt * derivative(t, y)
        y2 <- dt * derivative(t + dt/2, y + 0.5 * y1)
        y3 <- dt * derivative(t + dt/2, y + 0.5 * y2)
        y4 <- dt * derivative(t + dt,   y + y3)
        dy <- (y1 + 2 * y2 + 2 * y3 + y4)/6
        return (y + dy)
     }


FullY <- function(y,t,derivative,dt)

#----------------------------------------------------------------
# One runge-kutta step RK4, performed at time t
# Returns all intermediate steps
#----------------------------------------------------------------
     {
        # four steps
        y1 <- dt * derivative(t, y)
        y2 <- dt * derivative(t + dt/2, y + 0.5 * y1)
        y3 <- dt * derivative(t + dt/2, y + 0.5 * y2)
        y4 <- dt * derivative(t + dt,   y + y3)
        dy <- (y1 + 2 * y2 + 2 * y3 + y4)/6
        c(y,y+0.5*y1,y+0.5*y2,y+y3,y1/dt,y2/dt,y3/dt,y4/dt,y+dy)
     }

derivative <- function (x,y) mu*y
y   <- yini <- 0.01
Y   <- y
mu  <- 0.1
for (i in 1:50) {y<-NewY(y,i, derivative,1);Y<-c(Y,y)}

derplot <- function(x,y,dy,lab=NULL,py=0.05,px=0)
{
DD <- 2
points(x,y,cex=2,pch=16)
segments(x-DD,y-DD*dy,x+DD,y+DD*dy,lwd=2)
if(! is.null(lab)) text(x+px,y+py,lab)
}

par(mar=c(3,4,3,2))
x1  <- 30
y1  <- yini*exp(mu*x1)
dt  <- 10
YY  <- FullY(y1,1,derivative,dt)
curve(yini*exp(mu*x),25,45,lwd=1,xlab="",ylab="",ylim=c(0.,0.6),axes=FALSE,frame.plot=TRUE,main="Runge-Kutta method")
segments(30,-5,30,y1,lty=2)
segments(40,-5,40,YY[9],lty=2)
segments(35,0.05,35,yini*exp(mu*35),lty=2)

derplot(x1,YY[1],YY[5],"1")
derplot(x1+0.5*dt,YY[2],YY[6],"2",-0.02,0.7)
derplot(x1+0.5*dt,YY[3],YY[7],"3")
derplot(x1+dt,YY[4],YY[8],"4")

points(x1+dt,YY[9],cex=2,pch=16,col="darkgrey")
text(x1+dt+0.2,YY[9]-0.02,"estimated",col="darkgrey",adj=0)

segments (30,0.0,40,0.0)
text(35,0.02,expression(Delta~t))
segments (30,0.05,35,0.05)
text(32.5,0.07,expression(Delta~t/2))
axis(1,at=c(30,40),labels=c("t1","t2"))
writelabel("B")





# semi-implicit Euler
mu <- 0.9
fun <- function(x) exp(mu*x)
der <- function(x) mu*exp(mu*x)

dx  <- 1
x1  <- 1
x2  <- x1+dx

y1  <- fun(x1)
y2  <- fun(x2)

yse  <- y1 + dx*0.5*(der(x1) +der(x2))          # semi-implicit
yim <- y1 + dx* der(x2)                         # implicit

y1  <- fun(x1)
dx1 <- der(x1)
dx2 <- der(x2)

par(mar=c(3,4,3,2))
curve(fun(x),0,2.8,lwd=2,axes=FALSE,frame.plot=TRUE,xlab="",ylab="",
       ylim=c(0,8),main="Semi-Implicit euler method")
segments(x1,0,x1,y1,lty=2)

segments(x1,y1,x2,yse,lty=1)
segments(x2,0,x2,yse,lty=2)
dd <- 0.5
segments(x2-dd,y2-dd*dx2,x2+dd,y2+dd*dx2,lty=1)
segments(x1-dd,y1-dd*dx1,x1+dd,y1+dd*dx1,lty=1)
axis(2,at=c(-1,100),labels=NA,font=2,cex.axis=0.8)
axis(1,at=c(x1,x2),labels=c("t1","t2"))

points(x1,y1,pch=16,col="black",cex=2)
points(x2,yse,pch=16,col="darkgrey",cex=2)
points(x2,y2,pch=18,col="black",cex=2)
writelabel("D")
subtitle()

#####################################################
# FIG. 6.6. numerical approximation in space             #
#####################################################
par(las=0)
par(mfrow=c(2,2))
par(mar=c(2,2,2,1))
emptyplot(xlim = c(0,0.8),ylim = c(-0.25,0.8))
box(col="grey")
basemid <-  c(0.4,0)
dy <- 0.115
ds <- 0.05
straightarrow(c(0.4,0.75),c(0.4,-0.175),arr.pos=1,arr.type="triangle",arr.len=0.3)
for (i in 1:6)
filledcylinder(ry=0.2,rx=ds,len=dy, mid=basemid+c(0.,(i-1)*dy),angle=90,dr=0.01,
             col=grey(0.3+i*0.1),lwd=1,topcol=grey(0.3+i*0.1),lcol="black")
lines(c(0.4,0.4),c(0.75,0.625),lwd=2)
lines(c(0.4,0.4),c(0.625,-0.15),lwd=1,lty=2)
ii <- -0
text(0.72,-ii,expression(x[6]),cex=1.5)
text(0.72,-ii+dy,expression(x[5]),cex=1.5)
text(0.72,-ii+2*dy,expression(x[4]),cex=1.5)
text(0.72,-ii+3*dy,expression(x[3]),cex=1.5)
text(0.72,-ii+4*dy,expression(x[2]),cex=1.5)
text(0.72,-ii+5*dy,expression(x[1]),cex=1.5)
writelabel("A",at=-0.05,line=0)

emptyplot()
box(col="grey")
basemid <-  c(0.4,0.3)
dy <- 0.3
TXT <- c(expression(Delta~x[i+1]),expression(Delta~x[i]))
for (i in 1:2)
{
filledrectangle(wx=0.5,wy=dy,mid=basemid+c(0,(i-1)*dy),col="grey",lcol="black")
text(basemid[1]-0.05,basemid[2]+(i-1)*dy,TXT[i],cex=1.5,adj=0)
straightarrow(basemid+c(-0.1,(i-1.5)*dy),basemid+c(-0.1,(i-0.5)*dy),arr.pos=1,
              arr.type="triangle",arr.adj=1)
straightarrow(basemid+c(-0.1,(i-0.5)*dy),basemid+c(-0.1,(i-1.5)*dy),arr.pos=1,
              arr.type="triangle",arr.adj=1)
}

straightarrow(basemid+c(0.35,(0)*dy),basemid+c(0.35,(1.)*dy),arr.pos=1,
              arr.type="triangle",arr.adj=1)
straightarrow(basemid+c(0.35,(1)*dy),basemid+c(0.35,(0)*dy),arr.pos=1,
              arr.type="triangle",arr.adj=1)
text(basemid[1]+0.4,basemid[2]+(0.5)*dy,expression(Delta~x["i,i+1"]),cex=1.5,adj=0)

text(0.5,0.9,"Distances",cex=1.5)
writelabel("B",at=-0.05,line=0)


emptyplot()
box(col="grey")
basemid <-  c(0.5,0.3)
dy <- 0.2
txt <- c(expression(C[i+1]),expression(C[i+1]),expression(C[i]))
for (i in 1:2)
filledrectangle(wx=0.5,wy=dy,mid=basemid+c(0,(i-1)*dy),col="grey",lcol="black")
text(basemid[1]+0.325,basemid[2]+dy,"i",cex=1.5)
text(basemid[1]+0.325,basemid[2]   ,"i+1",cex=1.5)
for (i in 1:2) text(basemid[1],basemid[2]+(i-1)*dy,txt[i+1],cex=1.5)

text(0.5,0.9,"Concentrations",cex=1.5)
writelabel("C",at=-0.05,line=0)

emptyplot()
box(col="grey")
basemid <-  c(0.5,0.2)
dy <- 0.2
for (i in 1:3)
filledrectangle(wx=0.5,wy=dy,mid=basemid+c(0,(i-1)*dy),col="grey",lcol="black")
straightarrow(basemid+c(-0.1,2*dy-0.065),basemid+c(-0.1,dy+0.065),arr.pos=1,arr.type="triangle",arr.len=0.3)
straightarrow(basemid+c(-0.1,dy-0.065),basemid+c(-0.1,0.065),arr.pos=1,arr.type="triangle",arr.len=0.3)

straightarrow(basemid+c(0.1,2*dy-0.065),basemid+c(0.1,dy+0.065),arr.pos=1,arr.type="triangle",arr.len=0.3)
straightarrow(basemid+c(0.1,dy-0.065),basemid+c(0.1,0.065),arr.pos=1,arr.type="triangle",arr.len=0.3)
straightarrow(basemid+c(0.1,dy+0.065),basemid+c(0.1,2*dy-0.065),arr.pos=1,arr.type="triangle",arr.len=0.3)
straightarrow(basemid+c(0.1,0.065),basemid+c(0.1,dy-0.065),arr.pos=1,arr.type="triangle",arr.len=0.3)
#straightarrow(basemid-c(0,0.065),basemid-c(0,dy-0.065),arr.pos=1,arr.type="triangle",arr.len=0.3)
text(basemid[1]+0.325,basemid[2]+2*dy,"i-1",cex=1.5)
text(basemid[1]+0.325,basemid[2]+dy   ,"i",cex=1.5)
text(basemid[1]+0.325,basemid[2]      ,"i+1",cex=1.5)

text(basemid[1]-0.4,basemid[2]+1.5*dy,expression(Flux["i-1,i"]),cex=1.5)
text(basemid[1]-0.4,basemid[2]+0.5*dy,expression(Flux["i,i+1"]),cex=1.5)


text(0.5,0.9,"Fluxes",cex=1.5)
writelabel("D",at=-0.05,line=0)

subtitle()
####################################
# Fig. 6.7. Numerical diffusion
####################################

par(mfrow=c(1,1))
par(mar=c(5.1,4.1,4.1,2.1))
#----------------------#
# the model equations: #
#----------------------#

Growth<-function(t,Dens,parameters)
 {
     Num        <- c(0,Dens)
     Flux       <- growth*Num
     dDens      <- -diff(Flux)/delx
     list(dDens)     # result

  }  # end of model

#-----------------------#
# the model parameters: #
#-----------------------#
# dC/dt = -u*dC/dx
# C(0) = 10

growth   <- 1    # mm/day, the growth rate
delx     <- 1    # mm      thickness of boxes
numboxes <- 100
size     <- seq(0.5,100,times.out=numboxes)

# Initial condition
Numini   <- 2
Dens     <- c(rep(10,Numini),rep(0,(numboxes-Numini)))
State    <- c(Dens=Dens)

#------------------------------------#
# RUNNING the model: 100 size boxes  #
#------------------------------------#

times     <-seq(0,100,by=10)   # output wanted at these time intervals
out       <-ode(State,times,Growth,parms=0)

# the data in 'out' consist of: 1st col times,other columns: density

dens      <- out[,2:(numboxes+1)]
iout      <- c(1,seq(10,100,20))
#matplot(t(dens),type="l")

x0 <- c(0,9.999,10,11,11.0001,30)
td <- c(0,0,10,10,0,0)

dens10    <- dens [2,]


plot(size,dens10,type="l",ylim=c(0,10),xlim=c(0,20),lwd=2,ylab="density")
lines(x0,td,lwd=2,lty=2)

legend("topright",c("true size","modeled size"),lty=c(2,1),lwd=2)
subtitle()

################################################
# Fig. 6.8. Fiadeiro
################################################

Fiadeiro <- function (Pe) 
{
    sigma <- (1 + (1/tanh(Pe) - 1/Pe))/2
    return(sigma)
}
par(mfrow=c(1,1))
curve(Fiadeiro,0.01,100,log="x",lwd=2,ylab=expression(sigma),xlab="Pe",main="Fiadeiro weights")
subtitle()

################################################
# Fig. 6.9. Numerical errors in spatial discretisation
################################################
# general parameters
times     <- c(1/365,10/365,30/365,1/4,1/2,1)          # y    time at which output is required

# General parameters for all runs
DD        <- 0.1      # cm2/y   diffusion coefficient representing bioturbatin
u         <- 1       # cm/y    advection = sedimentation rate
init      <- 1000    # n.cm-3  initial density in upper layers
outtime   <- 6       # choose the output time point (as found in 'times') you want to see in the graphs
                     # default outtime=6 will show distribution after 1.0 year

# Function to calculate derivatives
Lumin <-function(t,LUMIN,parameters)
 {  with (as.list(c(LUMIN,parameters)),
    {
    FluxDiff     <- -DD*diff(c(LUMIN[1],LUMIN,LUMIN[numboxes]))/delx
    FluxAdv      <- u * (sigma * c(0,LUMIN) + (1-sigma) * c(0,LUMIN[2:numboxes],LUMIN[numboxes]))
    Flux         <- FluxDiff + FluxAdv

    dLUMIN       <- -diff(Flux)/delx
    
    list(dLUMIN )
     })
 }
 
# Function to run the model and produce output

Modelfunc <- function(ltt,diffm,numboxes)
{ 
    # complete set of parameters, depending on choices 
     delx      <- 3/numboxes                                     # thickness of boxes in cm
    if (diffm=="back")   sigma<-1
    if (diffm=="cent")   sigma<-0.5
    if (diffm=="fiad") {
                        if (DD>0) { Pe <- u * delx/DD
                                    sigma <- (1 + (1/tanh(Pe) - 1/Pe))/2
                                  }
                        if (DD==0) sigma <- 1
                       }
     pars<-as.list(c(          
              numboxes = numboxes,
              delx     = delx,
              sigma    = sigma,
              DD       = DD,
              u        = u
              ))

    # initial conditions
     borders         <- seq(from=0,by=delx,length.out=numboxes) 
     Ll              <- rep(0,times=numboxes)       # ind/m2
     for (i in 1:numboxes-1){     
       if (borders[i+1]<=0.5)                       Ll[i]  <- init
       if ((borders[i]<0.5) && (borders[i+1]>0.5))  Ll[i]  <- init * (0.5-borders[i])/delx
                            }

    # call
     state           <- c(LUMIN=Ll)
     out             <- ode.1D(state,times,Lumin,parms=pars,nspec=1)
    # store output
     Dl              <- out[outtime,2:(numboxes+1)]
    # output
     Distance  <- seq(from=delx/2, by=delx, length.out=numboxes) # distance of box centres from x=0
     if (diffm=="back") ttext<-"Backward Differences"
     if (diffm=="cent") ttext<-"Centered Differences"
     if (diffm=="fiad") ttext<-"Fiadeiro Scheme"
     if (ltt==1){ 
        plot(Dl,Distance,xlim=c(0,max(Dl)),ylim=c(numboxes*delx,0),type="l",
           main=ttext,xlab="Tracer density",ylab="Depth (cm)")
        legend(200,2, legend=c("n=300","n=120","n=60","n=30","n=15"),  
              lty=c(1,2,3,4,5),title="no. of boxes",pch=c(-1,-1,-1,-1,-1))
                } 
     if (ltt>1)lines(Dl,Distance,lty=ltt)
}

# Run model with number of different options
    par (mfrow=c(1,3),oma=c(1,1,4,1),mar=c(5, 4, 4, 2)+0.1,mgp=c(3,1,0)) 
    numboxlist <- c(300,120,60,30,15)
    
    for (diffm in c("back","cent","fiad")) {
     for (ltt in 1:5)                    {
       numboxes<-numboxlist[ltt]
       Modelfunc(ltt,diffm,numboxes)
                                         }
                                           }  
ttext<-paste("Luminophore distribution at t = ",round(times[outtime],2),"year")
mtext(outer=TRUE,side=3,ttext,cex=2)
subtitle()

###############################################################################
####======================================================================#####
####                            R case studies                            #####
####======================================================================#####
###############################################################################

## 6.10 Enzymatic reaction model ##

#-----------------------#
# the model parameters: #
#-----------------------#
                                     
parameters<-c(k1=0.01/24,            # parameter values
              k2=0.1/24,
              k3=0.1/24)

#-------------------------#
# the initial conditions: #
#-------------------------#
 
state     <-c(D=100,               
              I=10,
              E=1,
              F=1,
              G=0)

#----------------------#
# the model equations: #
#----------------------#

model<-function(t,state,parameters){
 with(as.list(c(state,parameters)),{  

    dD <- -k1*E*D + k2*I
    dI <-  k1*E*D - k2*I - k3*I*F
    dE <- -k1*E*D + k2*I + k3*I*F
    dF <-                - k3*I*F
    dG <-                  k3*I*F
    list(c(dD,dI,dE,dF,dG))          
    })
}


#----------------------#
# RUNNING the model:   #
#----------------------#

times     <-seq(0,300,0.5)         

out <-as.data.frame(ode(state,times,model,parameters))

par(mfrow=c(2,2), oma=c(0,0,3,0))   # set number of plots (mfrow) and margin size (oma)

plot (times,out$D,type="l",main="[D]",xlab="time, hours",ylab="mol/m3")
plot (times,out$F,type="l",main="[F]",xlab="time, hours",ylab="mol/m3")
plot (times,out$E,type="l",main="[E]",xlab="time, hours",ylab="mol/m3")
plot (times,out$I,type="l",main="[I]",xlab="time, hours",ylab="mol/m3")
mtext(outer=TRUE,side=3,"enzymatic reaction",cex=1.5)
subtitle()
 
################################################
# Fig. 6.11. Daphnia schema
################################################

par(mfrow=c(1,1))
par(mar=c(1,1,1,1),cex=0.8)

################################################
# drawing of a TRANSFER CULTURE
################################################
emptyplot()
dx <- 0.01
dy <- 0.1
xx <- c(0.0,1)
yy <- c(0.05,0.70) 
col <- grey(0.9)  
clock <- c(0.1,0.88)
clockrad<- 0.1    
lines(c(xx[1],xx[1],xx[2],xx[2]),c(yy[2],yy[1],yy[1],yy[2]))
rect(xx[1]+dx,yy[1]+dx,xx[2]-dx,yy[2]-dy,lwd=1,col=col,border=NA)

bentarrow (clock,c(xx[1]+0.25,xx[2]-0.25),arr.pos=1,arr.type="triangle",arr.length=0.5)
plotcircle(mid=clock,r=clockrad,lwd=1,col="white")
segments(clock[1],clock[2],clock[1],clock[2]+clockrad*0.8)
segments(clock[1],clock[2],clock[1]+clockrad*0.65,clock[2]-clockrad*0.3)

text(clock[1]+clockrad*1.65,clock[2]+0.025,"Food")

straightarrow(c(0.15,0.3),c(0.6,0.3),arr.pos=0)
textrect(mid=c(0.15,0.3),radx=0.075,lab="FOOD")
straightarrow(c(0.6,0.3),c(0.8,0.4))
straightarrow(c(0.6,0.3),c(0.8,0.2))
straightarrow(c(0.3,0.3),c(0.4,0.23),lwd=2,arr.pos=1)
straightarrow(c(0.4,0.3),c(0.5,0.37),lwd=2,arr.pos=1)
text(0.45,0.2,"Faeces")
text(0.55,0.4,expression(CO[2]))
bentarrow(from=c(0.8,0.4),to=c(0.935,0.3),arr.pos=1)
text(0.935,0.265,expression(CO[2]))
textrect(mid=c(0.8,0.4),radx=0.075,lab="BIOMASS")
textrect(mid=c(0.8,0.2),radx=0.075,lab="EGGS")
subtitle()

################################################
# Fig. 6.12. Daphnia model
################################################
## Daphnia model            ##

#----------------------#
# the model equations: #
#----------------------#


model<-function(t,state,parameters)
 {
 with(as.list(c(state)),{  # unpack the state variables

  # ingestion, size-dependent and food limited
  WeightFactor <- (IngestWeight-INDWEIGHT)/(IngestWeight-neonateWeight) 
  MaxIngestion <- maxIngest*WeightFactor      # /day     
  Ingestion    <- MaxIngestion*INDWEIGHT*FOOD / (FOOD + ksFood)

  Respiration  <- respirationRate * INDWEIGHT         # mugC/day 
  Growth       <- Ingestion*assimilEff - Respiration

  # Fraction of assimilate allocated to reproduction

  if (Growth <= 0. | INDWEIGHT<reproductiveWeight) Reproduction <- 0. 
  else {               # Fraction of growth allocated to reproduction.
    WeightRatio  <- reproductiveWeight/INDWEIGHT
    Reproduction <- maxReproduction * (1. - WeightRatio^2)
       }

  # rate of change 
  dINDWEIGHT <- (1. -Reproduction) * Growth
  dEGGWEIGHT <-      Reproduction  * Growth
  dFOOD      <- -Ingestion * numberIndividuals      

  # the output, packed as a list
    list(c(dINDWEIGHT, dEGGWEIGHT, dFOOD),   # the rate of change
         c(Ingestion    = Ingestion,            # the ordinary output variables
           Respiration  = Respiration,
           Reproduction = Reproduction))
    })

  }  # end of model

#----------------------#
# Moulting weight loss #
#----------------------#

Moulting   <- function ()

  {
   with(as.list(c(state)),{  # unpack the state variables
 
   # Relationship moulting loss and length
    refLoss   <-  0.24   #mugC
    cLoss     <-  3.1    #-

    # Weight lost during molts depends allometrically on the organism length
    INDLength    <- (INDWEIGHT /3.0)^(1/2.6)

    WeightLoss <- refLoss * INDLength^cLoss
    return(INDWEIGHT - WeightLoss)   # New weight
    })
  }

#-----------------------#
# the model parameters: #
#-----------------------#

neonateWeight      <-  1.1    #mugC
reproductiveWeight <-  7.5    #mugC
maximumWeight      <- 60.0    #mugC

ksFood             <- 85.0    #mugC/l
IngestWeight       <-132.0    #mugC
maxIngest          <-  1.05   #/day
assimilEff         <-  0.8    #-

maxReproduction    <-  0.8    #-
respirationRate    <-  0.25   #/day

# Dilution parameters !
transferTime       <-    2    # Days
foodInMedium       <-  509    # mugC/l

instarDuration     <-  3.0    # days
numberIndividuals  <-   32    #   -

#-------------------------#
# the initial conditions: #
#-------------------------#
 
state     <-c(
  INDWEIGHT = neonateWeight      , # mugC
  EGGWEIGHT = 0                  , # mugC    ! Total egg mass in a stage
  FOOD      = foodInMedium         # mugC
             )

#----------------------#
# RUNNING the model:   #
#----------------------#

TimeFrom     <- 0
TimeEnd      <- 40                          # duration of simulation, days
TimeMoult    <- TimeFrom + instarDuration   # next time (days) at which moulting 
TimeTransfer <- TimeFrom + transferTime     # next time (days) at which individuals are transferred

Time         <- TimeFrom
Outdt        <- 0.1                         # output time step
out          <- NULL                        # output array

while (Time < TimeEnd)
{
  TimeOut <- min(TimeMoult,TimeTransfer,TimeEnd)  # integrator runs till TimeOut
  times   <- seq(Time,TimeOut,by=Outdt)           # sequence of output times
  if (length(times)>1) {
  out1    <-as.data.frame(ode(state,times,model,parms=0))  # integrate
  out     <- rbind(out,out1)                                 # add output to output array
  lout    <- nrow(out1)                                      # last element of output
  state   <-c(
             INDWEIGHT = out1[lout,"INDWEIGHT"],
             EGGWEIGHT = out1[lout,"EGGWEIGHT"],
             FOOD      = out1[lout,"FOOD"])
  }
  if (Time >= TimeMoult)     # Moulting...
    {
     state[1]     <- Moulting()  # New weight individuals
     state[2]     <- 0.          # Put eggs = 0
     TimeMoult    <- Time +instarDuration          # next time at which moulting occurs
    }
  if (Time >= TimeTransfer)  # New medium...
    {
     state[3]     <- foodInMedium
     TimeTransfer <- Time + transferTime           # next time at which individuals are transferred
    }
  
  # Reset time, state variables
  Time   <- TimeOut 

 }

par(mfrow=c(2,2), oma=c(0,0,3,0))   # set number of plots (mfrow) and margin size (oma)
par(mar=c(5.1,4.1,4.1,2.1))
plot (out$time,out$FOOD        ,type="l",main="Food"              ,xlab="time, days",ylab="gC/m3")
plot (out$time,out$INDWEIGHT   ,type="l",main="individual weight" ,xlab="time, days",ylab="mugC")
plot (out$time,out$EGGWEIGHT   ,type="l",main="egg weight"        ,xlab="time, days",ylab="mugC")
plot (out$time,out$Ingestion   ,type="l",main="Ingestion"             ,xlab="time, days",ylab="mugC/day")

mtext(outer=TRUE,side=3,"DAPHNIA model",cex=1.5)
subtitle()

################################################
# Fig. 6.13. Daphnia relationships
################################################

par (mfrow=c(2,2))
curve(maxIngest*(IngestWeight-x)/(IngestWeight-neonateWeight),0,60,
      main="Max. ingestion rate",ylab="/d",xlab="ind. weight, mugC",lwd=2) 
curve(pmax(0., maxReproduction * (1. - (reproductiveWeight/x)^2)),0,60,
      main="fraction assimilate to reproduction ",ylab="-",
      xlab="ind. weight, muC",lwd=2) 
curve(((x /3.0))^(1/2.6),0,60,
      main="Individual length",ylab="mum",xlab="ind. weight, muC",lwd=2) 
curve(0.24*((x /3.0)^(1/2.6))^3.1,0,60,
      main="Weight loss during moulting",ylab="mug",xlab="ind. weight, mugC",lwd=2) 
subtitle()

################################################
# Fig. 6.14. 0-D zooplankton model
################################################

# Time and measured value of zooplankton concentration at sea boundary
# mg DWT/m3
fZooTime = c(0, 30,60,90,120,150,180,210,240,270,300,340,367)
fZooConc = c(20,25,30,70,150,110, 30, 60, 50, 30, 10, 20, 20)

# Traditional implementation
mod <- function (t,ZOO,parms,g,k)
{
  ZOOsea <- approx(fZooTime,fZooConc, xout=t)$y
   dZOO  <- k*(ZOOsea - ZOO) - g*ZOO
   return(list(dZOO,ZOOsea))
}

print(system.time(
Out <-ode(y=(ZOO=5),times= 0:365,func=mod,parms=NULL,k=0.015,g=0.05)
))


# The model integrated with Euler
euler <-function(start,end,delt,g,k)
 {
  times     <- seq(start,end,delt)
  nt        <- length(times)

  ZOOsea    <- approx(fZooTime,fZooConc, xout=times)$y

  ZOO       <- 5
  out       <- matrix(ncol=3,nrow=nt)

  for (i in 1:nt)
   {
     decay   <- g*ZOO
     input   <- k*(ZOOsea[i] - ZOO)

     out[i,] <- c(times[i],ZOO, input)

     dZOO    <- input - decay

     ZOO     <- ZOO + dZOO *delt
}

colnames(out) <- c("time","ZOO","input")
return(as.data.frame(out))

}
print(system.time(out <-euler(0,365,0.1,k=0.015,g=0.05)))

par (mfrow=c(1,1))
par (oma=c(0,0,0,2))

plot(fZooTime, fZooConc, type="b",xlab="daynr",ylab="mgDWT/m3",
     pch=15,main="Zooplankton model",lwd=2,ylim=c(0,150))
lines(out$time,out$ZOO,lwd=2,col="darkgrey")

legend("topright",c("Marine zooplankton","Estuarine zooplankton"),
       pch=c(15,NA),lwd=2,col=c("black","grey"))
subtitle()       


################################################
# Fig. 6.15. Aphid model
################################################
 

#----------------------#
# the model equations: #
#----------------------#
  
model <-function(t,APHIDS,parameters)
 {
    deltax     <- c (0.5,rep(1,numboxes-1),0.5)    
    Flux       <- -D*diff(c(0,APHIDS,0))/deltax
    dAPHIDS    <- -diff(Flux)/delx  + APHIDS*r

    # the output
      list(dAPHIDS )
  }  # end of model

#-----------------------#
# the model parameters: #
#-----------------------#

D         <- 0.3    # m2/day  diffusion rate
r         <- 0.01   # /day    net growth rate
delx      <- 1      # m       thickness of boxes
numboxes  <- 60 

Distance        <- seq(from=0.5,by=delx,length.out=numboxes)  # sequence, 1 m intervals

#--------------------------#
# Initial conditions:      #
#--------------------------#

APHIDS          <- rep(0,times=numboxes)       # ind/m2   The aphid density
APHIDS[30:31]   <- 1
state           <- c(APHIDS=APHIDS)            # the state variables are initialised
                  
#----------------------#
# RUNNING the model:   #
#----------------------#

times     <-seq(0,200,by=1)   # output wanted at these time intervals           
out       <- ode(state,times,model,parms=0)  # ode is integration routine

DENSITY   <- out[,2:(numboxes  +1)]

#------------------------#
# PLOTTING model output: #
#------------------------#

par(mfrow=c(1,1))
par(oma=c(0,0,3,0))   # set outer margin size (oma)
color= topo.colors

filled.contour(x=times,y=Distance,DENSITY,color= color,
               xlab="time, days", ylab= "Distance on plant, m",main="Density")
mtext(outer=TRUE,side=3,"Aphid model",cex=1.5)  # margin text
subtitle()


################################################
# Fig. 6.16. Aphid model
################################################
 

# 2. plot initial, intermediate and final densities, density versus time..
par(mfrow=c(2,2),oma=c(0,0,3,0))   # multiple figures on a row (2 rows, 2 cols), change margin size (oma)

plot(Distance,DENSITY[1,]  ,type="l",lwd=2,xlab="Distance, m",ylab="Density", main="initial condition")
plot(Distance,DENSITY[100,],type="l",lwd=2,xlab="Distance, m",ylab="Density", main="100 days")
plot(Distance,DENSITY[200,],type="l",lwd=2,xlab="Distance, m",ylab="Density", main="200 days")

meanAphid <- rowMeans(out[,2:ncol(out)])
plot(times,meanAphid  ,type="l",xlab="time, days",ylab="/m2",lwd=2, main="Density versus time") 

mtext(outer=TRUE,side=3,"Aphid model",cex=1.5)
subtitle()

################################################
# Fig. 6.17. Estuarine zooplankton
################################################

Zootran <-function(t,Zoo,pars)
{

  with (as.list(pars),{

    Flow   <- meanFlow+ampFlow*sin(2*pi*t/365+phaseFlow)
    seaZoo <- approx(fZooTime, fZooConc, xout=t)$y
    Input  <- +Flow * c(riverZoo, Zoo) +
              -Estar* diff(c(riverZoo, Zoo, seaZoo))
    dZoo   <- -diff(Input)/Volume + g *Zoo
    list(dZoo)
                       })
} 

#-----------------------#
# the forcing function: #
#-----------------------#

# Time and measured value of zooplankton concentration at sea boundary
# mg DWT/m3
fZooTime = c(0, 30,60,90,120,150,180,210,240,270,300,340,367)
fZooConc = c(20,25,30,70,150,110, 30, 60, 50, 30, 10, 20, 20)


#-----------------------#
# the model parameters: #
#-----------------------#

# the model parameters:  
pars   <-   c(riverZoo  = 0.0,          # river zooplankton conc
              g         =-0.05,         # /day  growth rate
              meanFlow  = 100*3600*24,  # m3/d, mean river flow
              ampFlow   = 50*3600*24,   # m3/d, amplitude
              phaseFlow = 1.4)          # -     phase of river flow

#--------------------------#
# Initialising morphology: #
#--------------------------#

# parameters defining the morphology
# cross sectional surface area is a sigmoid function of estuarine distance
nbox    <- 100                          
Length  <- 100000                           # m 

dx      <- Length/nbox                      # m

IntDist <- seq(0,by=dx,length.out=nbox+1)   # m
Dist    <- seq(dx/2,by=dx,length.out=nbox)  # m

IntArea <- 4000 + 76000 * IntDist^5 /(IntDist^5+50000^5)   # m2
Area    <- 4000 + 76000 * Dist^5    /(Dist^5+50000^5)      # m2

Volume  <- Area*dx                          # m3


#--------------------------#
# Transport coefficients:  #
#--------------------------#
# parameters defining the dispersion coefficients
# a linear function of estuarine distance

Eriver   <- 0                               # m2/d 
Esea     <- 350*3600*24                     # m2/d 
E        <- Eriver + IntDist/Length * Esea  # m2/d 

Estar  <- E * IntArea/dx                   # m3/d

#----------------------#
# RUNNING the model:   #
#----------------------#
ZOOP  <- rep(5,times=nbox)
times <- 1:365
out   <- ode.band(times=times,y=ZOOP,func=Zootran,parms=pars,nspec=1)

#------------------------#
# PLOTTING model output: #
#------------------------#

# Plot zooplankton; first get the data
par(mfrow=c(1,1))   # set margin size (oma)
par(oma=c(0,0,3,0))       # set margin size

color <- terrain.colors

filled.contour(x=times,y=Dist/1000,z=out[,-1],
               color= color,xlab="time, days",
               ylab= "Distance, km",main="Zooplankton, mg/m3")
mtext(outer=TRUE,side=3,"Marine Zooplankton in the Scheldt",cex=1.5) 
subtitle()

###############################################################################
####======================================================================#####
####                              R projects                              #####
####======================================================================#####
###############################################################################
################################################
# Fig. 6.18. Autocatalyse
################################################

# dilution plot

par(mfrow=c(1,2))
mar <- par(mar=c(2,2,2,2))
dd<-dilution(int=expression(d[r]))
text(dd$p2[1],dd$p2[2]+0.02,"Ain",font=2)
text(dd$p2[1],dd$p2[2]-0.02,"Bin",font=2)
elpos<-matrix(nr=5,data=c(0.45,0.65,0.37,0.55,0.73,  0.5,0.5,0.2,0.2,0.2))
tt<-treearrow(from=elpos[1:2,],to=elpos[3:5,],arr.side=2)
lab<-c("A","B","B","B","C")
text(0.55,0.378,"k")
for ( i in 1:5) textrect (elpos[i,],0.05,0.05,lab=lab[i],cex=1.5)
writelabel("A",line=-1,at=-0.05)


# model 
autocatalysis <- function(t,state,pars)
{

with (as.list(c(state,pars)),
{
dA <- dr*(Ain-A)-k*A*B
dB <- dr*(Bin-B)+k*A*B
dC <- -dr*C +k*A*B
return (list(c(dA,dB,dC)))
})

}
# model application

times <- seq(0,300,1)
state <- c(A=1,B=1,C=0)
parms <- c(Ain =1,Bin = 0.1, k = 0.05, dr = 0.05)
require(deSolve)
out   <- as.data.frame(ode(state,times,autocatalysis,parms))
ylim  <- c(0,1.1)

par(mar=mar)
plot(out$time,out$A,xlab="time",ylab="concentration",
      lwd=2,type="l",ylim=ylim,main="autocatalysis")
lines(out$time,out$B,lwd=2,lty=2)
lines(out$time,out$C,lwd=2,lty=3)

legend("topright",c("[A]","[B]","[C]"),lwd=2,lty=c(1,2,3)) 
subtitle()
################################################
# Fig. 6.19. Algae-DIN in dilution
################################################

par(mfrow=c(1,1))
par(mar=c(0,0,0,0))
dd<-dilution(int="dilrate")
text(dd$p2[1],dd$p2[2],"Nin",font=2)
elpos<-matrix(nr=2,data=c(0.4,0.7,0.3,0.3))
tt<-straightarrow(from=elpos[1,],to=elpos[2,])
lab<-c("DIN","Algae")
for ( i in 1:2) textrect (elpos[i,],0.08,0.05,lab=lab[i],cex=1.25)
 

subtitle()

par(ask=opar$ask)
par(mar=opar$mar)
par(oma=opar$oma)
par(mfrow=opar$mfrow)

