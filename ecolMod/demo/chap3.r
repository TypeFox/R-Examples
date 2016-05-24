##############################
## Soetaert and Herman      ##
## ecological modelling     ##
## Figures from chapter 3   ##
## Space and transport      ##
##############################

opar <- par()
par(mfrow=c(1,1))

par(ask=TRUE)
figNr <- 1
subtitle <- function()
{
 mtext(side=1,outer=TRUE,"Soetaert and Herman - chapter 3  ",cex=0.7,adj=1,line=-1.5)
 mtext(side=1,outer=TRUE,paste("  Fig. 3.",figNr,sep=""),cex=0.7,adj=0,line=-1.5)
 figNr <<- figNr +1
}

###############################################################################
####======================================================================#####
####                                Theory                                #####
####======================================================================#####
###############################################################################

################################################
# Fig. 3.1. Microscopic, macroscopic model
################################################

par(mar=c(2,2,2,2))
emptyplot(xlim=c(-0.95,0.95),ylim=c(-1.05,1.25),asp=FALSE)
mid <- c(-0.6,0.6)
rx=0.1
ry=0.2
filledcylinder(rx=rx,ry=ry,angle=90,col="white",lcolint="grey",mid=mid,lcol="black")
x<-(runif(4000)-0.5)*0.035+mid[1]
y<-runif(4000)-0.6+mid[2]
ijmin <- 0.1
ii <- which(y>ijmin)
points(x[ii],y[ii],pch=".")
mtext(side=3,at=-0.6,"initial", adj=0.5,font=3)


#filledcylinder(rx=0.01,ry=0.02,angle=90,col="black",topcol="black",mid=mid,lcol="black")

# second cylinder, random points
mid <- c(0.6,0.6)
filledcylinder(rx=rx,ry=ry,angle=90,col="white",lcolint="grey",mid=mid,lcol="black")
x<-rnorm(10000,0,0.3)*0.2+mid[1]
y<-runif(10000)*1.1-0.6+mid[2]
ii <- which (x<mid[1]-0.2 | x>mid[1]+0.2)

# don't worry about the NaNs message...
ij <- which(y<mid[2]-0.5)
Xa <- acos((x[ij]-mid[1])/ry)
ijmin <-mid[2]-0.5-rx*sin(Xa)   # remove all those
ikmin <-mid[2]-0.5-0.2*rx*sin(Xa)   # and part of those
ik <- which(y[ij]<ijmin)
ii<-c(ii,ij[ik])
ik <- which(y[ij]<ikmin)
ii <- c(ii,ij[ik[1:(length(ik)*0.8)]])

points(x[-ii],y[-ii],pch=".")
# top ellipse
plotellipse(rx=ry,ry=rx,col="white",mid=mid+c(0,0.5))
mtext(side=3,at= 0.6,"final", adj=0.5,font=3)
# arrows and text
segments(-0.2,0.6,0.2,0.6,lty=2)
Arrows(-0.2,0.6,0.2,0.6,arr.type="triangle",segment=FALSE)
text(0.0,1.2,"microscopic", adj=0.5,font=2)
text(0.0,1.15,"model", adj=0.5,font=2)

text(0.0,0.65,"time", adj=0.5,font=3)

Arrows(-0.6,-0.1,-0.6,-0.3,arr.type="triangle",arr.col="darkgrey")
Arrows(0.6,-0.1,0.6,-0.3,arr.type="triangle",arr.col="darkgrey")
text(0.0,-0.2,"statistical  averaging", adj=0.5,font=3)

# macro initial
Arrows(-0.6,-0.95,-0.6,-0.45,arr.type="triangle",lty=2)
Arrows(-0.8,-0.95,-0.35,-0.95,arr.type="triangle",code=3)
text(-0.85,-0.7, "density",srt=90)
lines(c(-0.62,-0.62,-0.58,-0.58),c(-0.95,-0.6,-0.6,-0.95),lwd=2)

# macro final
Arrows(0.6,-0.95,0.6,-0.45,arr.type="triangle",lty=2)
Arrows(0.4,-0.95,0.85,-0.95,arr.type="triangle",code=3)
text(0.35,-0.7, "density",srt=90)

# normal dist
x <- seq(-1,1,length.out=1000)
y <- dnorm(x,0,0.3)
x <- 0.6 +x*0.2
y <- 0.1*y-0.95
lines(x,y,lwd=2)

segments(-0.2,-0.7,0.2,-0.7,lty=2)
Arrows(-0.2,-0.7,0.2,-0.7,arr.type="triangle",segment=FALSE)
text(0.0,-0.45,"macroscopic", adj=0.5,font=2)
text(0.0,-0.5,"model", adj=0.5,font=2)
text(0.0,-0.65,"time", adj=0.5,font=3)
# boxes
rect(-0.925,-0.025,0.925,1.235)
rect(-0.925,-0.4,0.925,-1.05)

subtitle()

################################################
# Fig. 3.2. Discrete spatial model
################################################

par(mfrow=c(2,2))
par(mar=c(3,3,3,3))
require(deldir)
# landscape models
dd <- trunc(dim(volcano)/2)
VV <- matrix(nr=dd[1],nc=dd[2],data=0)
for ( i in 1:dd[1])                                                             {
 for (j in 1:dd[2]) VV[i,j] <- 0.25*(
 volcano[(i-1)*2+1,(j-1)*2+1]+volcano[(i-1)*2+2,(j-1)*2+1]+
 volcano[(i-1)*2+1,(j-1)*2+2]+volcano[(i-1)*2+2,(j-1)*2+2])
                    }
 image( VV, col = greycol(100),axes=FALSE)
 grid(dd[1],dd[2],col="darkgrey",lty=1,lwd=0.5)
 box()
writelabel("A")

# patch model
emptyplot()
m1 <- c(0.2,0.8)
m2 <- c(0.6,0.6)
m3 <- c(0.4,0.2)
m4 <- c(0.8,0.15)
straightarrow(m1,m2)
straightarrow(m1,m3)
straightarrow(m3,m2)
straightarrow(m4,m3)
filledcircle(mid=m1, r1=0.1,col=grey(0.9))
filledellipse(mid=m2, rx1=0.1,ry1=0.15,angle=100,col=grey(0.8))
filledellipse(mid=m2, rx1=0.1,ry1=0.15,angle=200,col=grey(0.8))
filledellipse(mid=m3, rx1=0.05,ry1=0.15,angle=200,col=grey(0.95))
filledmultigonal(mid=m3, rx=0.12, col=grey(0.925),nr=8.5)
filledmultigonal(mid=m4, rx=0.08, col=grey(0.925),nr=5.2)

box()
writelabel("B")


# triangulation model
# the original plot.deldir from package deldir is strongly simplified...
plot.deldir_simple <- function (x)
{
    n  <- x$n.data
    plot(0, 0, type = "n", xlim = x$rw[1:2], ylim = x$rw[3:4], xlab = "",
            ylab = "",  axes = FALSE)
    segments(x$dirsgs[,1], x$dirsgs[,2], x$dirsgs[,3], x$dirsgs[,4],lty = 2)
    points(x$summary[1:n, 1], x$summary[1:n, 2], pch = 1)
}

 try <- deldir(x=runif(100),y=runif(100),list(ndx=2,ndy=2),c(0,1,0,1))
 plot.deldir_simple(try) 
 box()
writelabel("C")

# cellular automaton
par(mar=c(3,3,3,3))

radical <- function(mid)
{
 dy <- 0.05
 dx <- 0.04
 segments(mid[1]-dx,mid[2]+dy,mid[1],mid[2],lwd=4)
 segments(mid[1],mid[2],mid[1]+dx,mid[2]+dy*2,lwd=2)
}

particle <- function(mid) points(mid[1],mid[2],pch=16,cex=4)
cross    <- function(mid) points(mid[1],mid[2],pch="x",cex=3)

plot(0,type="n",xlim=c(0,1.4),ylim=c(0,1.4),axes=FALSE,xlab="",ylab="")
for (i in seq(0.2,1.2,0.2)) abline(h=i)
for (i in seq(0.2,1.2,0.2)) abline(v=i)
radical(c(0.3,1.25))
radical(c(1.1,1.25))
radical(c(0.3,0.45))
radical(c(0.5,0.25))
particle(c(0.3,1.1))
particle(c(0.3,0.3))
particle(c(0.5,0.5))
particle(c(0.9,0.5))
particle(c(1.1,0.3))
particle(c(0.9,1.1))
particle(c(1.1,1.1))
cross(c(1,1.1))
cross(c(1.1,0.5))
arrows(0.3,0.3,0.43,0.3,length=0.1,lwd=2)
arrows(0.5,0.5,0.35,0.5,length=0.1,lwd=2)
arrows(0.9,0.5,1.05,0.5,length=0.1,lwd=2)
arrows(0.9,1.1,1.05,1.1,length=0.1,lwd=2)
arrows(1.1,1.1,1.1,1.23,length=0.1,lwd=2)
arrows(0.3,1.1,0.3,1.23,length=0.1,lwd=2)
arrows(1.1,0.3,1.1,0.45,length=0.1,lwd=2)
rect(-0.1,-0.1,0.4,0.2,col="white")
text(-0.025,0.1,"\\sr",vfont=c("serif","plain"))
text(-0.01,0.1,"allowed",adj=0)
text(-0.025,0.0,"x")
text(-0.01,0.0,"not allowed",adj=0)
box()
writelabel("D")
subtitle()

################################################
# Fig. 3.3. zone-of-influence
################################################

par(mfrow=c(1,2))
par(mar=c(2,2,2,2))

# 1. radiusses

plot(c(0, 0), type = "n", ylab = "", asp = 1, xaxt = "n",
       yaxt = "n", frame.plot = FALSE, xlim = c(0,0.4),
       ylim = c(0.7,0.95),main="zone of influence",sub="",xlab="")
       
plotcircle(r=0.13,mid=c(0.14,0.85),lwd=1)
plotcircle(r=0.1 ,mid=c(0.3,0.85),lwd=1)

x1<-getellipse(rx=0.1,ry=0.1 ,mid=c(0.3,0.85),from=pi*0.7,to=pi*1.3)
x2<-getellipse(rx=0.13,ry=0.13,mid=c(0.14,0.85),from=-pi*0.21,to=pi*0.21)

polygon(x=c(x1[,1],x2[,1],x1[1,1]),y=c(x1[,2],x2[,2],x1[1,2]),col="darkgrey")
points(0.3,0.85,pch=16,cex=1.25)
points(0.14,0.85,pch=16,cex=1.25)

arrows(0.3,0.85,0.3+0.1,0.85,length=0.1,lwd=2)
arrows(0.14,0.85,0.14-0.13,0.85,length=0.1,lwd=2)
text(0.115,0.82,"zone of influence")
text(0.115,0.795,"radius")
segments(0.24,0.85,0.24,0.72)
text(0.23,0.71,"zone of overlap")

text(0.3,0.92,"B",cex=1.5)
text(0.14,0.92,"A",cex=1.5)

# overlap

par(mar=c(5,4,3,2))

dtheta1 <- pi*(1.8-1.2)
dtheta2 <- pi*(0.76-0.23)

r1     <- 0.1
area1 <- dtheta1/2/pi*pi*r1^2
tria1 <- cos(dtheta1/2)*r1*sin(dtheta1/2)*r1


r     <- 0.13
area2 <- dtheta2/2/pi*pi*r^2
tria2 <- cos(dtheta2/2)*r*sin(dtheta2/2)*r

totarea <- area1+area2 -tria1-tria2
overlapB <- totarea/pi/r1^2
overlapA <- totarea/pi/r^2

# relative performance
plot(x=c(0,0.9,1),y=c(1,0,0),lwd=2, type = "l", ylab = "-",
     xlim = c(0,1),main="relative performance",xlab="relative overlap")

points(overlapA,1-1/0.9*overlapA,pch=16,cex=1.5)
text(overlapA,1-0.9*overlapA+0.1,"A")
points(overlapB,1-1/0.9*overlapB,pch=16,cex=1.5)
text(overlapB,1-0.9*overlapB+0.1,"B")
subtitle()


################################################
# Fig. 3.6  Well-stirred tank, 0-D transport
################################################
nf <- layout(matrix(c(1,1,2,3),2,2, byrow=TRUE), c(2,2), c(4,2))

#par(mfrow=c(2,2))
par(mar=c(3,3,3,3))
dd<-dilution(main=c("Stock","Stirred tank"),int="Flow,Q")
text(dd$p2[1],dd$p2[2],expression(A["in"]),font=2,cex=1.2)
text(dd$p1[1],dd$p1[2]+0.03,"Volume V",font=2)
text(dd$p1[1],dd$p1[2]-0.03,"[A]",font=2)

box()
writelabel("A",at=-0.05)
par(mar=c(1,3,1,1))
col<-grey(0.9)
x  <- seq(-0.9*pi,0.9*pi,len=100)
x2 <- seq(-0.8*pi,0.8*pi,len=100)
plot(x,1-cos(x),axes=FALSE,type="n",xlab="",ylab="",xlim=range(x)*1.5,ylim=c(-0.5,3.0))
rect(-4,1.6,4,1.9,col=col,border=NA)
lines(c(-4,4),c(1.6,1.6),lwd=1)
lines(c(-4,4),c(1.9,1.9),lwd=1)
polygon(x2,1-cos(x2),col=col,border=NA)
xx<-seq(-0.71*pi,0.71*pi,len=100)
lines(xx,1-cos(xx),lwd=1)
Arrows(-3.5,1.75,-2.5,1.75,arr.type="triangle",arr.length=0.35)
Arrows(2.5,1.75,3.5,1.75,arr.type="triangle",arr.length=0.35)

selfarrow(c(-0.3,1),curve=0.25,arr.type="triangle",arr.length=0.35,lwd=1)
selfarrow(c(0.3,1),curve=0.25,arr.pos=0.001,path="R",
         arr.type="triangle",arr.length=0.35,lwd=1)

box()
writelabel("B")

par(mar=c(1,1,1,3))

plot(0,axes=FALSE,type="n",xlab="",ylab="",xlim=c(-1,1),ylim=c(-1,1.0))
rect(-0.8,-0.8,0.8,0.6,col=col,border=NA)
lines(c(-0.8,0.8),c(-0.8,-0.8))
lines(c(-0.8,-0.8),c(-0.8,0.8))
lines(c(0.8,0.8),c(-0.8,0.8))

Arrows(-0.1,0.45,-0.1,0.75,arr.type="triangle",arr.adj=1)
Arrows(0.1,0.75,0.1,0.45,arr.type="triangle",arr.adj=1)

selfarrow(c(0.0,-0.1),curve=0.2,arr.pos=0.001,path="H",
         arr.type="triangle",arr.length=0.35,lwd=1)

box()
writelabel("C")

subtitle()

################################################
# Fig. 3.5  General transport
################################################

par(mar=c(1,1,1,1))
plot(0,type="n",xlim=c(-1,1),ylim=c(-0.6,0.5),axes=FALSE,xlab="",ylab="")
col <- grey(seq(0.2,1,length.out=100))
col<-c(col,rev(col))
cex<-1.75

filledcylinder (rx=0.15,ry=0.4,len=1,col=col,lcol="black",lwd=1,lcolint=grey(0.25),
                lwdint=1,ltyint=3,topcol=grey(0.5),delt=1.15)

segments(-1,0,-0.5,0)

segments(0.5,0,1,0)
Arrows(-0.8,0,-0.5,0,arr.type="triangle",arr.length=0.5,lwd=5,arr.adj=1)
Arrows(0.5,0,0.8,0,arr.type="triangle",arr.length=0.5,lwd=3,arr.adj=1)

text(0.0,0.5,expression(Delta~V),cex=cex*0.9)
text(-0.5,0.225,expression(A[x]),cex=cex)
text(0.5,0.225,expression(A[x+Delta~x]),cex=cex)

text(-0.75,0.065,expression(J[x]),cex=cex)
text(0.85,0.065,expression(J[x+Delta~x]),cex=cex)

segments(-0.5,0,-0.5,-0.5,lty=3,col=grey(0.25))
segments(0.5,0,0.5,-0.5,lty=3,col=grey(0.25))

text(-0.5,-0.55,expression(x),cex=cex)
text(0.5,-0.55,expression(x+Delta~x),cex=cex)
subtitle()

################################################
# Fig. 3.6. Advection and diffusion graph
################################################

par(mar=c(1,4,3,1))#,mfrow=c(2,2))
nf<-layout(matrix(nr=2,nc=2,byrow=TRUE,c(1,2,3,3)))
#layout.show(nf)

openplotmat()

my <-0.9       #position y
rx <- 0.075    #"long
ry <- 0.125
len <- 0.6
dd <- 0.6
col <-grey(0.95)

cols <- greycol(100,interval=c(0.05,0.7))
plotellipse(mid=c(0.5-0.5*len,my),from=pi,to=2*pi,rx=rx,ry=ry,lcol=cols[100],lwd=2)
L <- len/100
for (i in 1:100)
 cylindersegment(mid=c(0.2+(i-1)*L,my),rx=rx,ry=ry+ry*(i-1)*dd/100,to=3/2*pi,from=pi,col=cols[101-i],len=L,delt=1+dd/100,lwd=2)
plotellipse(mid=c(0.5+len/2,my),from=pi,to=2*pi,rx=rx*(1+dd),ry=ry*(1+dd),col=cols[1],lwd=1)

xx <- seq(0.15,0.65,by=0.1)
x2 <- xx+0.045
Arrows(xx,0.885,x2,0.885,arr.length=0.2,lwd=2,arr.type="triangle")

filledrectangle(mid=c(0.5,0.3),wx=0.5,wy=0.4,col=cols,lcol="grey")
yy <- matrix(nc=2,data=runif(100))
yy[,1] <- 0.28 + yy[,1]*0.45
yy[,2] <- pmax(0.2,0.4 + (0.5-rexp(50,2))*0.15)
#seq(0.15,0.5,by=0.08)
dy <- yy-0.045
points(yy)
Arrows(yy[,1],yy[,2],yy[,1],yy[,2]-0.045,arr.length=0.1,lwd=1,arr.type="triangle")
box()

writelabel("A")

# Diffusion graph

xy <- matrix(nc=2,data=runif(100))
xy[,1]<-xy[,1]+0.1
xy[,2]<-xy[,2]-0.1

dxy <- xy + matrix(nc=2,data=(runif(100)-0.5)*0.2)
plot(xy,type="n",axes=FALSE,xlab="",ylab="",xlim=c(0,2.9),ylim=c(-1.8,1))
Arrows(xy[,1],xy[,2],dxy[,1],dxy[,2],arr.len=0.175,arr.type="triangle")
#rect(-0.1,-0.1,1.1,1.1)
rect(0,-0.2,1.2,1.0)
#text(0.,1.025,"A",cex=1.75)

xy <- matrix(nc=2,data=runif(200),byrow=TRUE)+c(rep(1.5,100),rep(0,100))
xy[,1]<-xy[,1]+0.1
xy[,2]<-xy[,2]-0.1

points(xy,cex=0.75)
pos <- matrix(nc=2,byrow=TRUE,c(1.75,0.75,2.5,0.5,2.5,0.1,1.8,0.3,2.8,0.9,2.2,0.75))*0.9
pos[,1]<-pos[,1]+0.1
pos[,2]<-pos[,2]-0.1

for (i in 1:nrow(pos)) plotellipse(mid=pos[i,],rx=rx<-max(0.1,0.2*runif(1)),ry=rx,arr.len=0.25,
angle=runif(1)*360,arr.type="triangle",arrow=TRUE,arr.pos=1,from=0,to=3*pi/2)
#rect(1.3,-0.1,2.8,1.1)
rect(1.4,-0.2,2.9,1.)

rect(0,-1.8,2.9,-0.4)
#text(0,-0.475,"C",cex=1.75)

xl = 0.3
xr = 2.6

yd = -1.5
yu = -0.7

segments(xl,yd,xr,yd,lwd=3)
segments(xl,yu,xr,yu,lwd=3)

xl = 0.5
xr = 1.5
yu = -0.8
yd = -1.4

xx <- c(0.7,0.825,0.9, 0.95, 0.9,0.825,0.7)+0.5
yy<- yd +(0:7)/10#+0.05
for ( i in 1:7)
 Arrows(xl,yy[i],xx[i],yy[i],arr.type="triangle",arr.len=0.2)
 
plotellipse(mid=c(1.55,-1.35),rx=0.05,ry=0.05,arr.len=0.15,lwd=1,
           angle=180,arr.type="triangle",arrow=TRUE,arr.pos=1,from=0,to=3*pi/2)

plotellipse(mid=c(1.55,-0.85),rx=0.05,ry=0.05,arr.len=0.15,lwd=1,
           angle=-180,arr.type="triangle",arrow=TRUE,arr.pos=0,from=0,to=3*pi/2)
#box(col="grey")
box()
writelabel("B")

# Diffusion/advection of a dye spill - analytical model solutions 

Ct <- function(t,x,E,v,k) C0*exp(k*t)*exp(-(((x-L)-v*t)^2)/4/E/t)/2/sqrt(pi*E*t)

xx<-seq(1,100,len=1000)
v <-1     #velocity
E <-0.25  #dispersion
k <- -0.0 #decay
C0 <- 5   #conc of initial spill
L  <- 20  #pos of initial spill
tt <- 20  #time 

#initial condition
Cin<-rep(0,100);Cin[L]<-5

# advection +diffusion 
CC <-Ct(tt,xx,E,v,k)

# advection only 
C10<-rep(0,100);C10[tt*v+L]<-5*exp(k*tt)

# diffusion only
CE <-Ct(tt,xx,E,v=0,k)

#windows(5,5)
par(mar=c(3,7,3,1))
openplotmat(ylim=c(0.1,0.9))
col<- greycol(100,c(0,1))
filledrectangle(mid=c(0.4,0.8),wx=0.1,wy=0.5,angle=270,col=intpalette(col,x.to=Cin/max(Cin)),lcol="black")
filledrectangle(mid=c(0.4,0.6),wx=0.1,wy=0.5,angle=270,col=intpalette(col,x.to=C10/max(C10)),lcol="black")
filledrectangle(mid=c(0.4,0.4),wx=0.1,wy=0.5,angle=270,col=intpalette(col,x.to=CE/max(CE)),lcol="black")
filledrectangle(mid=c(0.4,0.2),wx=0.1,wy=0.5,angle=270,col=intpalette(col,x.to=CC/max(CC)),lcol="black")

Arrows(0.25,0.6,0.31,0.6,arr.col="darkgrey",lcol="darkgrey")
Arrows(0.25,0.2,0.31,0.2,arr.col="darkgrey",lcol="darkgrey")
curvedarrow(c(0.22,0.47),c(0.28,0.47),curve=-0.2,lwd=1,arr.pos=0.85,arr.length=0.25,arr.adj=0,arr.col="darkgrey",lcol="darkgrey")
curvedarrow(c(0.28,0.47),c(0.22,0.47),curve=0.2,lwd=1,arr.pos=0.85,arr.length=0.25,arr.adj=0,arr.col="darkgrey",lcol="darkgrey")
curvedarrow(c(0.32,0.27),c(0.38,0.27),curve=-0.2,lwd=1,arr.pos=0.85,arr.length=0.25,arr.adj=0,arr.col="darkgrey",lcol="darkgrey")
curvedarrow(c(0.38,0.27),c(0.32,0.27),curve=0.2,lwd=1,arr.pos=0.85,arr.length=0.25,arr.adj=0,arr.col="darkgrey",lcol="darkgrey")

text(0.05,0.8,"t=0",adj=0)
text(0.05,0.6,"t=20",adj=0)
text(0.05,0.4,"t=20",adj=0)
text(0.05,0.2,"t=20",adj=0)

text(0.75,0.8,"dye spill",adj=0)
text(0.75,0.6,"advection only",adj=0)
text(0.75,0.4,"diffusion only",adj=0)
text(0.75,0.22,"advection +",adj=0)
text(0.75,0.18,"diffusion",adj=c(0,0.5))
box()
writelabel("C")
subtitle()

################################################
# Fig. 3.7. diffusion/advection of estuary/river/lake
################################################

par(mar=c(1,1,2,1))
par(mfrow=c(2,1))
openplotmat(main="river,estuary")
col<- greycol(100,c(0,1))
my <-0.7       #position y
rx <- 0.075    #"long
ry <- 0.25
len <- 0.6
dd <- 0.6
plotellipse(mid=c(0.5-0.5*len,my),from=pi,to=2*pi,rx=rx,ry=ry,col="grey",lwd=2)

cylindersegment(mid=c(0.5,my),rx=rx,ry=ry,to=3/2*pi,from=pi,col="grey",len=len,delt=1+dd,lwd=2)
for (i in seq(0,1,len=9)) plotellipse(mid=c(0.5-0.5*len+len*i,my),from=pi,to=1.5*pi,rx=rx*(1+i*dd),ry=ry*(1+i*dd),lcol="darkgrey",lwd=1)
plotellipse(mid=c(0.5-0.5*len,my),from=pi,to=1.5*pi,rx=rx,ry=ry,lcol="black",lwd=2)
plotellipse(mid=c(0.5+len/2,my),from=pi,to=2*pi,rx=rx*(1+dd),ry=ry*(1+dd),col="darkgrey",lwd=2)

Arrows(0.3,0.8,0.7,0.8)
box(col="grey")
xx <- seq(-0.2*pi,1.2*pi,len=20)
yy <- sin(xx)

#par(mar=c(3,5,3,5))
x  <- seq(-0.9*pi,0.9*pi,len=100)
x2 <- seq(-0.8*pi,0.8*pi,len=100)
plot(x,1-cos(x),axes=FALSE,type="n",xlab="",ylab="",main="lake",xlim=range(x)*1.2,ylim=c(-0.1,2.0))
polygon(x2,1-cos(x2),col="grey")
for (y in seq(0,1.7,by=0.25)) lines(c(-acos(1-y),acos(1-y)),c(y,y),col="darkgrey")
lines(x,1-cos(x),lwd=2)

Arrows(2.8,1.5,2.8,0.2)
box(col="grey")
subtitle()

################################################
# Fig. 3.8. 1-D shapes 
################################################

par(mar=c(2,2,2,2),mfrow=c(2,2))

# FIRST the CONSTANT SURFACE MODEL
emptyplot(c(-1,1),frame.plot=TRUE,main="parallel isosurfaces")
filledcylinder(rx=0.15/2,ry=0.4/2,angle=90,col=grey(0.6),lcol="black",
               lcolint=grey(0.55),mid=c(0.5,0),botcol=grey(0.55))
for ( i in seq(-0.5,0.5,by=0.05)) plotellipse  (rx=0.2,ry=0.15/2,mid=c(0.5,i),col=grey(0.6),lwd=1)
plotellipse  (rx=0.2,ry=0.15/2,mid=c(0.5,0.5),col="white",lwd=2)
Arrows(0.5,-0.6,0.5,-0.7,arr.type="triangle",arr.length=0.2)
segments(0.5,0.7,0.5,-0.5,col=grey(0.5))
filledcylinder(rx=0.15/2,ry=0.2,angle=90,col=grey(0.6),lcol="black",
               lcolint=grey(0.55),mid=c(-0.5,0),botcol=grey(0.55))
Arrows(-0.5,0.7,-0.5,-0.7,arr.type="triangle",arr.length=0.2)
segments(-0.5,0.5,-0.5,-0.5,col=grey(0.5))
text(-0.5,0.3,adj=0,"M",font=2,cex=1.5)
text(-0.5,0.15,adj=0,"o",font=2,cex=1.5)
text(-0.5,0.0,adj=0,"d",font=2,cex=1.5)
text(-0.5,-0.15,adj=0,"e",font=2,cex=1.5)
text(-0.5,-0.3,adj=0,"l",font=2,cex=1.5)
text(-0.5,-0.8,"x",cex=1.5)
text(-0.4,0.65,"constant surface",cex=1.2,adj=0)
segments(-0.45,0.55,-0.4,0.6)
writelabel("A",at=0.1,line=-2)

# CYLINDRICAL MODEL 1
plot(0,type="n",xlim=c(-0.9,0.9),ylim=c(-1,1),axes=FALSE,xlab="",ylab="",
        frame.plot=TRUE,main="cylindrical isosurfaces")
col <- grey(seq(0.2,1,length.out=100))
filledcylinder(rx=0.15/2,ry=0.35/2,angle=90,col=grey(0.6),lcol="black",
               lcolint=grey(0.55),mid=c(0.5,0),botcol=grey(0.55))
for ( i in seq(0.2,0.8,by=0.2)) plotellipse  (rx=0.175*i,ry=0.15/2*i,mid=c(0.5,-0.5),lcol=grey(0.45),lwd=1)
for ( i in seq(0.2,0.8,by=0.2)) plotellipse  (rx=0.175*i,ry=0.15/2*i,mid=c(0.5,0.5),lcol=grey(0.45),lwd=1)
for ( i in seq(0.2,0.8,by=0.2)) segments  (0.5-0.175*i,-0.5,0.5-0.175*i,0.5,col=grey(0.45),lwd=1)
for ( i in seq(0.2,0.8,by=0.2)) segments  (0.5+0.175*i,-0.5,0.5+0.175*i,0.5,col=grey(0.45),lwd=1)
Arrows(0.5,0.5,0.8,0.5,arr.type="triangle",arr.length=0.2)
filledcylinder(rx=0.15/2,ry=0.35/2,angle=90,col=grey(0.6),lcol="black",
               lcolint=grey(0.55),mid=c(-0.5,0),botcol=grey(0.55))
Arrows(-0.5,0.5,-0.2,0.5,arr.type="triangle",arr.length=0.2)
text(-0.15,0.5,"r",cex=1.5,adj=0)
text(-0.5,0.3,adj=0,"M",font=2,cex=1.5)
text(-0.5,0.15,adj=0,"o",font=2,cex=1.5)
text(-0.5,0.0,adj=0,"d",font=2,cex=1.5)
text(-0.5,-0.15,adj=0,"e",font=2,cex=1.5)
text(-0.5,-0.3,adj=0,"l",font=2,cex=1.5)
writelabel("B",at=0.1,line=-2)

# CYLINDRICAL MODEL 2: cylindrical coordinates on a flat surface
plot(0,type="n",xlim=c(-1.2,1.2),ylim=c(-1.2,1.2),axes=FALSE,xlab="",ylab="",
asp=TRUE,frame.plot=TRUE,main="circular isolines")
plotcircle (0.7,col=grey(0.6),mid=c(0,0))
for ( i in seq(0.05,0.65,by=0.05)) plotcircle(i,lwd=1,lcol=grey(0.2))
Arrows(0.,0.,0.8,0.,arr.type="triangle",arr.length=0.2)
text(0.9,0.0,"r",cex=1.5,adj=0)
writelabel("C",at=0.1,line=-2)


## Spherical coordinates
emptyplot(c(-1.75,1.75),frame.plot=TRUE,main="spherical isosurfaces")
rseq <- seq(0.4,1,length.out=100)
col  <- grey(rev(rseq))
filledellipse (rx1=1,ry1=1.0,col=col)
ry <- 0.35
rss<-1
pow <- 1
plotellipse(rx=1,ry=ry,angle=-90,from=pi,to=2*pi,col=grey(0.95),lwd=1)
plotellipse(rx=1,ry=0.0,angle=90,from=pi,to=2*pi,lwd=1,lcol="grey")
plotellipse(rx=1,ry=rss,angle=90,from=pi,to=2*pi,col=grey(0.95),lwd=1)
for (i in seq(0.1,0.9,by=0.1))plotellipse(rx=i,ry=ry*i^pow,angle=-90,from=pi,to=2*pi,lwd=0.5)
for (i in seq(0.1,0.9,by=0.1))plotellipse(rx=i,ry=rss*i,angle=90,from=pi,to=2*pi,lwd=1)
Arrows(0,0,1.2,0,arr.type="triangle",arr.length=0.2)
text(1.2,0.2,"r",cex=1.5)
writelabel("D",at=0.1,line=-2)
subtitle()

################################################
# Fig. 3.9  porous media
################################################

par(mar=c(0,0,0,0),mfrow=c(1,1))
emptyplot(c(0,1),c(0,1.2),asp=FALSE)

rect(0.5,0.2,1.0,1.1,col=grey(0.95))
ngrain <- 5000
pos <- matrix(ncol=2,data=runif(2*ngrain))
expd <- rexp(ngrain,1.5)
ii   <- which(expd<1)
pos[ii,2] <-expd[ii]
pos[,1] <- 0.52+pos[,1]*0.46
pos[,2] <- 0.22+pos[,2]*0.66
angle <- 360*runif(ngrain)
dcol <- 0.1
bcol <- 0.4
siz  <- 0.012#5
for (i in 1:ngrain) filledrectangle(mid=pos[i,],wx=siz,wy=siz,angle=angle[i],
                    col=grey(bcol+runif(1)*dcol))
Arrows(code=3,0.1,0.2,0.1,0.875,arr.adj=1,arr.type="triangle",arr.length=.1)
Arrows(code=3,0.1,0.885,0.1,1.1,arr.adj=1,arr.type="triangle",arr.length=.1)

text(0.35,1,expression(phi==1))
text(0.35,0.55,expression(phi<1))

text(0.075,0.55,"sediment",srt=90,font=3)
text(0.075,1.0,"water",srt=90,font=3)

rect(0.10,0.05,0.30,0.15,col=grey(0.95))
    for (i in 1:100) filledrectangle(mid=c(0.11,0.06)+runif(2)*c(0.18,0.08),wx=0.015,wy=0.015,angle=angle[i],
                    col=grey(bcol+runif(1)*dcol))

rect(0.45,0.05,0.65,0.15,col=grey(0.95))
rect(0.80,0.05,1.0,0.15,col=grey(bcol))
    for (i in 1:100) filledrectangle(mid=c(0.81,0.06)+runif(2)*c(0.18,0.08),wx=0.015,wy=0.015,angle=angle[i],
                    col=grey(bcol+runif(1)*dcol))

text(0.2,0.02,"bulk sediment")
text(0.55,0.02,"liquid phase")
text(0.9,0.02,"solid phase")

dd <- seq(0,4,0.1)
x <- c(1,1,0.72+0.28*exp(-1*dd))-0.55
y <- c(1.1,0.9,seq(0.9,0.2,length.out=length(x)-2))
lines(x,y,lwd=2)
rect(0.15,0.2,0.48,1.1)
text((0.15+0.48)*0.5,1.15,"Porosity")
subtitle()


################################################
# Fig. 3.10. 1-D shapes - boundaries
################################################

par(mar=c(2,2,2,2),mfrow=c(2,2))

# FIRST the CONSTANT SURFACE MODEL
emptyplot(c(-0.6,0.6),c(-1,1),frame.plot=TRUE,asp=FALSE,main="parallel isosurfaces")
filledcylinder(rx=0.15/2,ry=0.4/2,angle=90,col="white",topcol=grey(0.6),lcol="black",
               lcolint=grey(0.55),mid=c(-0.2,0),botcol=grey(0.55),len=1.2)
segments(-0.15,-0.6, 0.3,-0.1)
segments(-0.15, 0.6, 0.3, 0.1)
Arrows(-0.2,0.8,-0.2,-0.8,arr.type="triangle",lty=2,arr.length=0.2)
text(-0.2,-0.9,"x",cex=1.5)
text(0.25,0,"boundaries",cex=1.2,adj=0)
writelabel("A",at=0.1,line=-2)

# Cylindrical coordinates
plot(0,type="n",xlim=c(-0.6,0.6),ylim=c(-1,1),axes=FALSE,xlab="",ylab="",
        frame.plot=TRUE,main="cylindrical isosurfaces")
filledcylinder(rx=0.15/2,ry=0.4/2,angle=90,col=grey(0.8),lcol="black",
               lcolint=grey(0.55),mid=c(-0.1,0),botcol=grey(0.55),len=1.2)
Arrows(-0.1,0.6,-0.4,0.6,arr.type="triangle",lty=2,arr.length=0.2)
text(-0.5,0.6,"r",cex=1.5,adj=0)
segments(-0.1,0.6,-0.1,-0.6,lwd=2)
segments(-0.1,-0.3,0.3,-0.1)
segments(0,0.3,0.3,0.1)
text(0.25,0,"boundaries",cex=1.2,adj=0)
writelabel("B",at=0.1,line=-2)

# cylindrical coordinates on a flat surface
plot(0,type="n",xlim=c(-1.2,1.2),ylim=c(-1.2,1.2),axes=FALSE,xlab="",ylab="",
asp=TRUE,frame.plot=TRUE,main="cylindrical isolines")
plotcircle (0.5,col="white",mid=c(-0.1,0))
Arrows(-0.1,0.,-0.85,0,arr.type="triangle",lty=2,arr.length=0.2)
points(-0.1,0,pch=18,bg="black")

segments(-0.1,-0.5,0.4,-0.5)
segments(-0.1,0,0.4,-0.3)
text(0.4,-0.4,"boundaries",cex=1.2,adj=0)

text(-1,0.0,"r",cex=1.5,adj=0)
writelabel("C",at=0.1,line=-2)


## Spherical coordinates
emptyplot(c(-1.75,1.75),frame.plot=TRUE,main="spherical isosurfaces")
rseq <- seq(0.4,1,length.out=100)
col  <- grey(rev(rseq))
filledellipse (rx1=1,ry1=1.0,col=col)
ry <- 0.35
rss<-1
pow <- 1
plotellipse(rx=1,ry=ry,angle=-90,from=pi,to=2*pi,col=grey(0.95),lwd=1)
plotellipse(rx=1,ry=rss,angle=90,from=pi,to=2*pi,col=grey(0.95),lwd=1)
plotellipse(rx=1,ry=0.0,angle=90,from=pi,to=2*pi,lwd=1,lcol="darkgrey")
Arrows(0,0,1.2,0,arr.type="triangle",lty=2,arr.length=0.2)
text(1.2,0.2,"r",cex=1.5)
segments(-0.6,-0.5,0.5,-1.2)
segments(0,0,0.7,-1.1)
points(0,0,pch=18,col="black")
text(0.55,-1.3,"boundaries",cex=1.2,adj=0)

writelabel("D",at=0.1,line=-2)
subtitle()

################################################
# Fig. 3.11. Boundaries and donut
################################################

par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
NumLayers <- 100
color <-gray(seq(1,0.3,length.out=50))
color <- c(rev(color),color)
maxrad <- 1.0
emptyplot(xlim=c(-2,2),ylim=c(-2,2.5))
y1<-1.95
dy <- 0.4
dx<- 0.4
DD <- dx/dy
straightarrow(c(-1.9,y1+0.2),c(-1.9,y1-0.1),arr.type="triangle",arr.pos=1)
straightarrow(c(-0.4,y1+0.2),c(-0.4,y1-0.1),arr.type="triangle",arr.pos=1)
straightarrow(c(-1.1,y1),c(-1.1,y1-0.3),arr.type="triangle",arr.pos=1)
text(-1.9,y1-0.2,"x",font=2)
text(-0.4,y1-0.2,"x",font=2)
text(-1.1,y1-0.4,"x",font=2)

polygon(x=c(-2.1,-0.6,-0.6+dx,-2.1+dx),y=c(y1,y1,y1+dy,y1+dy),col="lightgray")
lines(x=c(-2.1,-0.6,-0.6+dx),y=c(y1,y1,y1+dy),col="black",lwd=2)

for (i in seq(0.1,0.4,by=0.1)) lines(c(-2.1+DD*i,-0.6+DD*i),c(y1+i,y1+i))
for (i in seq(-2.1,-0.6,by=0.1)) lines(c(i,i+dx),c(y1,y1+dy))

polygon(x=c(0.2,1.7,1.7+dx,0.2+dx),y=c(y1,y1,y1+dy,y1+dy),col="lightgray")
lines(x=c(0.2,1.7,1.7+dx),y=c(y1,y1,y1+dy),col="black",lwd=2)

for (i in seq(0.1,0.4,by=0.1)) lines(c(0.2+DD*i,1.7+DD*i),c(y1+i,y1+i))
for (i in seq(0.2,1.7,by=0.1)) lines(c(i,i+dx),c(y1,y1+dy))
plotellipse(0.09,0.07,mid=c(1.3,y1+0.325),arrow=TRUE,
arr.length=0.3,arr.type="triangle",from=-pi/3,to=pi*1.2,arr.pos=1)

plotellipse(0.09,0.07,mid=c(0.5,y1+0.22),arrow=TRUE, angle=60,
arr.length=0.3,arr.type="triangle",from=-pi/3,to=pi*1.2,arr.pos=1)

plotellipse(0.09,0.07,mid=c(1.,y1+0.07),arrow=TRUE, angle=180,
arr.length=0.3,arr.type="triangle",from=-pi/3,to=pi*1.2,arr.pos=1)

plotellipse(0.09,0.07,mid=c(1.8,y1+0.2),arrow=TRUE, angle=240,
arr.length=0.3,arr.type="triangle",from=-pi/3,to=pi*1.2,arr.pos=1)

#straightarrow(c(-1.1,2.4),c(-1.1,2.0),arr.type="triangle",arr.pos=1)

mid <- c(0,0.6)
plotellipse(rx=1.8,ry=0.33,mid=mid,from=0.77*pi ,to=2.225*pi,lwd=2,arrow=TRUE,
arr.pos=0.5,arr.len=0.4)

y<- c(mid[2]+0.1,mid[2]+0.4)
rect(-1.5,y[1],1.5,y[2],col="lightgray")
for (i in seq(y[1],y[2],by=0.1)) lines(c(-1.5,1.5),c(i,i))
for (i in seq(-1.5,1.5,by=0.1)) lines(c(i,i),y)

mid <- c(0,mean(y)-0.06)

plotellipse(rx=0.15,ry=0.3,mid=mid,from=pi/4,to=1.9*pi,lwd=2)
ell <- getellipse(rx=0.15,ry=0.3,mid=mid,from=pi*0.99,to=1.11*pi)
Arrows(ell[1,1],ell[1,2],ell[nrow(ell),1],ell[nrow(ell),2],lwd=2)

filledellipse(mid=c(0,-0.9),rx1=1,rx2=0.5,col=color)
rect(-2.2,-2.0,2.2,1.35,border="grey")
rect(-2.2,1.4,-0.05,2.5,border="grey")
rect(0.05,1.4,2.2,2.5,border="grey")
text(-1.8,1.2,"C",cex=1.5)
text(0.25,2.4,"B",cex=1.5)

text(-1.95,2.4,"A",cex=1.5)

subtitle()

################################################
# Fig. 3.12. 1-D biogeochemical model boundaries
################################################

par(mfrow=c(1,2))
par(mar=c(2,4,3,1))
# require shape

sedimentfig <- function(upbnd1,upbnd2,main)
{

plot(0,type="n",xlim=c(-1,0),ylim=c(-1,0.8),axes=FALSE,xlab="",ylab="",
     main=main)
tpos <- -0.6# -0.35
adj <- 0

model<- expression(frac(partialdiff~C,partialdiff~t)==-v*frac(partialdiff~C,partialdiff~x) +
            D*frac(partialdiff^2~C,partialdiff~x^2)-k*C)
filledcylinder(rx=0.1,ry=0.15,angle=90,col="white",lcol="black",
               lcolint="grey",mid=c(-0.8,0),topcol=grey(0.5)) 
#Arrows(-0.8,0.8,-0.8,0.5,arr.type="triangle",arr.length=0.25,lwd=3,arr.adj=1)               
segments(-0.8,0.5,-0.8,-0.5,col=grey(0.5),lty=2)
Arrows(-0.8,-0.5,-0.8,-0.8,arr.type="triangle",lty=2,arr.length=0.2)               

text(-0.8,-0.95,"x",cex=1.5)
text(tpos,0,model,adj=0,cex=1)
rect(tpos-0.025,-0.15,0,0.15)


text(tpos,0.68,upbnd1,cex=1.,adj=adj)
text(tpos,0.55,upbnd2,cex=1.2,adj=adj)
text(tpos,-0.68,"Concentration boundary",cex=1.,adj=adj)
text(tpos,-0.9,expression (group("",C,"|")[infinity] == 0),cex=1.2,adj=adj)
text(tpos,-0.9,expression (group("",C,"|")[infinity] == 0),cex=1.2,adj=adj)

}

sedimentfig("Flux boundary",expression (group("",J,"|")[0] == flux0),"particulate substance")
Arrows(-0.8,0.8,-0.8,0.5,arr.type="triangle",arr.length=0.8,lwd=5,arr.adj=1)  
box(col="grey")
writelabel("A")

sedimentfig("Concentration boundary",expression (group("",C,"|")[0] == C0),"dissolved substance")
box(col="grey")
writelabel("B")
subtitle()

################################################
# Fig. 3.13. File spatial_omexdia.r
################################################

N      <- 100
Depth  <- seq(0.05,by=0.1,len=100)
out    <- OMEXDIAsteady()

# Steady-state concentrations in sediment
CONC  <- out$steady

O2    <- CONC[(2*N+1):(3*N)]
NO3   <- CONC[(3*N+1):(4*N)]
NH3   <- CONC[(4*N+1):(5*N)]

ii <- min(which(O2<0.1))
opd <- Depth[ii]
ii <- min(which(NO3<0.001))
npd <- Depth[ii]
mxO2 <- max(O2)

par(mfrow=c(1,2),mar=c(4,2,2,2))
plot(O2,Depth,ylim=c(10,0),type="n",axes=FALSE,xlab="",ylab="Depth")
rect(0,0,mxO2,11,col=NA,border="grey")

rect(0,opd,mxO2,npd,col="lightgrey",border=NA)
rect(0,npd,mxO2,11,col="darkgrey",border=NA)

lines(O2,Depth,lwd=2,col="grey")
par(new=TRUE)
plot(NO3,Depth,ylim=c(10,0),type="l",lwd=2,col=grey(0.4),axes=FALSE,xlab="",ylab="")
par(new=TRUE)
plot(NH3,Depth,ylim=c(10,0),type="l",lwd=2,axes=FALSE,xlab="",ylab="")
legend ("bottom",bg="lightgrey",col=c("grey",grey(0.4),"black"),lwd=2,legend=c("Oxygen","Nitrate","Sulphide"))


plot(O2,Depth,ylim=c(10,0),type="n",axes=FALSE,xlab="",ylab="Depth")
rect(0,0,255,11,col=NA,border="grey")
rect(0,opd,255,npd,col="lightgrey",border=NA)
rect(0,npd,255,11,col="darkgrey",border=NA)

text(1,0.3,"Oxic zone",adj=0,cex=1.2)
text(10,0.8,"oxic mineralisation, nitrification",adj=0,cex=0.9)

text(1,2.8,"Suboxic zone",adj=0,cex=1.2)
text(10,3.3,"denitrification",adj=0,cex=0.9)

text(1,8,  "Anoxic zone",adj=0,cex=1.2)
text(10,8.5,"anoxic mineralisation",adj=0,cex=0.9)
subtitle()

################################################
# Fig. 3.14. 1-D boundaries 2-layered model
################################################

par(mar=c(2,4,3,1))

sedimentfig <- function(model1,model2,intbnd,main,ci=TRUE)
{
plot(0,type="n",xlim=c(-1.1,0),ylim=c(-1,0.8),axes=FALSE,xlab="",ylab="",
     main=main)
tpos <- -0.6# -0.35
adj <- 0
filledcylinder(rx=0.1,ry=0.15,angle=90,col="white",lcol="black",
               lcolint="grey",mid=c(-0.8,0),topcol=grey(0.5)) 
text(-0.9,0.25,"I",cex=1.2)               
text(-0.9,-0.25,"II",cex=1.2)               
plotellipse(rx=0.15,ry=0.1,mid=c(-0.8,0))
text(tpos,0.4,model1,adj=adj,cex=1)
text(tpos,-0.4,model2,adj=adj,cex=1)
rect(tpos-0.025,0.25,0,0.55)
rect(tpos-0.025,-0.25,0,-0.55)
Arrows(-0.8,0.8,-0.8,0.5,arr.type="triangle",arr.length=0.8,lwd=5,arr.adj=1)               
text(-0.775,0.7,"F0",adj=0,cex=1)
segments(-0.8,0.5,-0.8,-0.5,col=grey(0.5),lty=2)
Arrows(-0.8,-0.5,-0.8,-0.8,arr.type="triangle",lty=2,arr.length=0.2)               

text(tpos,0.05,expression (C[I] == C[II]),cex=1.2,adj=adj)
text(tpos,-0.05,intbnd,cex=1.2,adj=adj)

text(-0.8,-0.95,"x",cex=1.5)
text(tpos,0.68,expression (group("",J,"|")[0] == F0),cex=1.2,adj=adj)
if(ci) text(tpos,-0.9,expression (group("",C,"|")[infinity] == 0),cex=1.2,adj=adj)
box(col="grey")

}

model1<- expression(frac(partialdiff~C,partialdiff~t)==-v*frac(partialdiff~C,partialdiff~x) +
            Db*frac(partialdiff^2~C,partialdiff~x^2)-k*C)
model2<- expression(frac(partialdiff~C,partialdiff~t)==-v*frac(partialdiff~C,partialdiff~x) -k*C)
intbnd   <- expression (J[I] == J[II])
main <- "2-layered model of C-dynamics"
sedimentfig (model1,model2,intbnd,main,ci=FALSE)


writelabel("A")

model2<- model1
intbnd   <- expression (J[I]+FI == J[II])
main <- "Non-local exchange model"
sedimentfig (model1,model2,intbnd,main)
Arrows(-1.1,0.,-0.95,0.,arr.type="triangle",arr.length=0.8,lwd=5,arr.adj=1)               
segments(-1.1,0.8,-1.1,0.0,lwd=5)
text(-1.075,0.7,"FI",adj=0,cex=1)
writelabel("B")
subtitle()

################################################
# Fig. 3.15. Autocatalysis, dilution
################################################

par(mfrow=c(1,2))
par(mar=c(0,0,0,0))
dd<-dilution(int=expression(d[r]))
text(dd$p2[1],dd$p2[2]+0.02,"Ain",font=2)
text(dd$p2[1],dd$p2[2]-0.02,"Bin",font=2)
elpos<-matrix(nr=5,data=c(0.45,0.65,0.37,0.55,0.73,  0.5,0.5,0.2,0.2,0.2))
tt<-treearrow(from=elpos[1:2,],to=elpos[3:5,],arr.side=2)
lab<-c("A","B","B","B","C")
text(0.55,0.378,"k")
for ( i in 1:5) textrect (elpos[i,],0.05,0.05,lab=lab[i],cex=1.5)

# Autocatalysis model
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

out   <- as.data.frame(ode(state,times,autocatalysis,parms))
ylim  <- range(c(out$A,out$B,out$C))
par(mar=c(5.1,4.1,4.1,2.1))
plot(out$time,out$A,xlab="time",ylab="concentration",
      lwd=2,type="l",ylim=ylim,main="autocatalysis")
lines(out$time,out$B,lwd=2,lty=2)
lines(out$time,out$C,lwd=2,lty=3)

legend("topright",c("A","B","C"),lwd=2,lty=c(1,2,3)) 
subtitle()


###############################################################################
####======================================================================#####
####                            R case studies                            #####
####======================================================================#####
###############################################################################

################################
## Micro- macroscopic model   ##
## of random movement         ##
################################

par(mfcol=c(2,2))

#======================
# microscopic model
#======================

nind    <- 100 
nsteps  <- 100
pos     <- matrix (nrow=nsteps+1,ncol=nind,data=0)

# random steps, between -0.5 and 0.5
for (i in 1:nsteps) pos[i+1,] <- pos[i,]+runif(nind)-0.5


matplot(pos,type="l",col="black",lty=1,xlab="step",ylab="Position")
mtext("microscopic model",side=3,line=2,cex=1.2)

# density plot
posmicro <- matrix(ncol=nsteps+1,nrow=100)
for (i in 1:(nsteps+1)) 
  posmicro[,i] <- density(pos[i,],from=-10,to=10,n=100)$y


yy     <- 1:nsteps
xx     <- seq(-10,10,length.out=100)

matplot(posmicro[,seq(1,100,by=4)],type="l",lty=1,col="black",
        axes=FALSE,frame.plot=TRUE,ylab="Density")


#=============================
# macroscopic model
#=============================
# estimate the parameters of the continuous time model:
vari <- apply(X=pos,MARGIN=1,FUN=var)
l1   <- lm(vari~ c(1:101) +0)

Ds  <-  coef(l1)/2   #0.075/2   # diffusion coefficient
ini <- 1       # initial condition


# analytical solution of the 1-D diffusion equation in cartesian coordinates  
xx <-seq(-10,10,length=100)
tt <-seq(1,100,by=1)
posmacro <- outer(xx,tt,
    FUN =function (xx,tt) ini/(2*sqrt(pi*Ds*tt))*exp(-xx^2/(4*Ds*tt)))

persp(xx,tt,z=posmacro,theta=150,box=TRUE,axes=TRUE,border=NA,
     xlab="position",ylab="time",zlab="density",col=drapecol(posmacro))  
mtext("macroscopic model",side=3,line=2,cex=1.2)

matplot(posmacro[ ,seq(1,100,by=4)],type="l",lty=1,col="black",
        axes=FALSE,frame.plot=TRUE,ylab="Density")
subtitle()

################################
## DIFFUSION in cellular      ##
## automaton                  ##
################################

#=====================================================

Diffuse <- function (Particles)
#-----------------------------------------------------
# Performs one diffusion step
#-----------------------------------------------------

  {
  
  # select the positions where there are Particles 
  x      <- which(Particles>0)
  
  # if random number > 0.5, particle moves to right (+1) if < 0.5 to left (-1)
  rnd            <- runif(length(x))   
  move           <- c(x[rnd<0.5]-1,x[rnd>=0.5]+1)

  # wrap around edges if necessary 
  move[move<1]     <- move[move<1]    +ncell
  move[move>ncell] <- move[move>ncell]-ncell  

  # Mark the new positions that are not empty              
  notfree <- which (Particles[move]>0)

  # movement allowed only when there is just one element moves to new position             
  duplo  <- which(move%in%move[duplicated(move)])
  free   <- move[-c(duplo,notfree)]
  source <- c(x[rnd<0.5],x[rnd>=0.5])[-c(duplo,notfree)]

  # source particles are emptied, new positions become 1
  Particles[source] <- 0
  Particles[free]   <- 1  
  return(Particles) 

} # END SUBROUTINE Diffuse

#=====================================================



##############################
# Initialising               #
##############################
  
ncell      <- 200
Particles  <- rep (0,ncell) 
Particles [c(10:30,50:70,90:110,130:150,170:190)] <- 1
nsteps     <- 1000

##############################
# RUNNING the model:         #
##############################

# perform nsteps diffusion steps 

Grid <- matrix(ncol=nsteps, nrow=ncell)

for (j in 1:nsteps) {
  Particles <- Diffuse(Particles)
  Grid[,j]  <- Particles
         }   

##############################
# PLOTTING model output:     #
##############################


par(oma = c(0,0,3,0),mfrow=c(1,2))
image (y=1:nsteps,z=Grid, col = c(0,1), axes=FALSE, 
       ylim=c(nsteps,1),xlab="",ylab="")
box()
mtext(outer=TRUE,side=3,"Diffusion in cellular automaton",cex=1.5)
subtitle()   


################################
## Plant competition          ##
## Function embedded in DLL   ##
## generated with Fortran     ##
################################


# competition function in FORTRAN
competition <- function(cells,nstep=100)

{
 nspec   <- nrow(replacement)
 ncell   <- nrow(cells)
 cells   <- matrix(nrow=ncell,ncol=ncell,data=as.integer(cells))
 sumdens <- matrix(nrow=nstep,ncol=nspec,as.integer(0))
 seed    <- runif(1)*-10
 res <- .Fortran("lattice",nspec=nspec,ncell=ncell,nstep=as.integer(nstep),
                 cells=cells,replacement=as.double(replacement),
                 sumdens=sumdens,seed=as.integer(seed))
 
 return(list  (cells=matrix(nr=ncell,nc=ncell,res$cells),
               density=res$sumdens))

}


##############################################
# application
##############################################
species <- c("Lolium","Agrostis","Holcus","Poa","Cynosurus")

replacement <- matrix(ncol=5,byrow=TRUE,data=c(
1.0,0.02,0.06,0.05,0.03,
0.23,1.0,0.09,0.32,0.37,
0.06,0.08,1.0,0.16,0.09,
0.44,0.06,0.06,1.0,0.11,
0.03,0.02,0.03,0.05,1.0 ) )


ini   <-c(4,5,1,3,2)
cells <- matrix(40,40,data=0)

cells[,1:8]  <-ini[1] ; cells[,9:16] <-ini[2] ;
cells[,17:24]<-ini[3] ; cells[,25:32]<-ini[4]
cells[,33:40]<-ini[5]

nstep=100
A100  <- competition(cells,nstep=100)
A200  <- competition(A100$cells,nstep=100)


# graphs

par(mfrow=c(2,2),mar=c(2,2,2,2))
col   <-c("grey","lightblue","blue","darkblue","black")

image(cells,col=col,zlim=c(1,5),axes=FALSE,main="initial")
text(x=rep(0.1,5),y=seq(0.1,0.9,length.out=5),
     labels=species[ini],col="white",adj=0,font=2)
image(A100$cells,col=col, zlim=c(1,5),axes=FALSE,
      main="100 steps")
image(A200$cells,col=col, zlim=c(1,5),axes=FALSE,
      main="200 steps")
matplot(rbind(A100$density,A200$density),type="l",lwd=2,lty=1,
        col=col,xlab="time",ylab="",axes=FALSE,frame.plot=TRUE)

mtext(outer=TRUE,side=3,"Spatial competition model",cex=1.5)
subtitle()

par(ask=opar$ask)
par(mar=opar$mar)
par(oma=opar$oma)
par(mfrow=opar$mfrow)



