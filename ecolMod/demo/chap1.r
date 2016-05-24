
##############################
## Soetaert and Herman      ##
## ecological modelling     ##
## Figures from chapter 1   ##
## Introduction             ##
##############################

opar <- par()
par(ask=TRUE)

figNr <- 1
subtitle <- function()
{
 mtext(side=1,outer=TRUE,"Soetaert and Herman - chapter 1  ",cex=0.7,adj=1,line=-1.5)
 mtext(side=1,outer=TRUE,paste("  Fig. 1.",figNr,sep=""),cex=0.7,adj=0,line=-1.5)
 figNr <<- figNr +1
}


###############################################################################
####======================================================================#####
####                                Theory                                #####
####======================================================================#####
###############################################################################



# THIS example is toggled off - to avoid needing to install mapdata
par(mar=c(2,0,2,0))
par(mfrow=c(2,2))

#in the book it uses 'mapdata ' here just maps...
#require(mapdata)                                  
#m<-map('worldHires', c('Belgium','Netherlands')) 
#map.text('worldHires','Belgium',add=TRUE)   

require(maps)
m<-map('world', c('Belgium','Netherlands'))
map.text('world','Belgium',add=TRUE)   
text(5.5,52,"Netherlands")
map.axes()

box(col="grey")
writelabel("A",at=0.055,line=-1.55)

openplotmat()
polygon(c(0.35,0.25,0.65,0.55),c(0.7,0.3,0.3,0.7),col="lightgrey")
polygon(c(0.1,0.25,0.65,0.8),c(0.45,0.3,0.3,0.55),col="grey")
box(col="grey")
writelabel("B",at=0.055,line=-1.55)

names <- c("DETRITUS","MEIO","MACRO",expression(NH[3]))
M <- matrix(nrow=4,ncol=4,byrow=TRUE,data=c(
#   d m M  n
    0,"","", 0, #d
    1,0,0, 0,#m
    2,5,0, 0, #M
    3,4,6,0  #n
    ))

pp<-plotmat(M,pos=c(1,2,1),curve=0,name=names,lwd=1,my=0.0,cex.txt=0.7,
            box.lwd=2,box.size=0.15,box.type="square",box.prop=0.5,box.cex=0.9,
            arr.type="triangle",arr.pos=0.6,shadow.size=0.01,arr.len=0.25)

# extra arrows: flow 5 to Detritus and flow 2 to detritus
nh3     <-pp$comp[4,]
bentarrow(nh3+c(0.15,0),nh3+c(0.2,-0.1),arr.pos=1,arr.type="triangle",arr.len=0.25)

box(col="grey")
writelabel("C",at=0.055,line=-1.55)
subtitle()
  
########################################
# fig 1.2: Zooplankton energy balance
########################################

par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
openplotmat()
bentarrow(from=c(0.5,0.5),to=c(0.85,0.3),arr.pos=1,arr.type="triangle")
bentarrow(from=c(0.5,0.5),to=c(0.15,0.3),arr.pos=1,arr.type="triangle")
straightarrow(from=c(0.5,0.9), to=c(0.5,0.65),arr.pos=1,arr.adj=1,arr.type="triangle")
straightarrow(from=c(0.5,0.5), to=c(0.5,0.1),arr.pos=1,arr.adj=1,arr.type="triangle")
textrect(mid=c(0.5,0.5),radx=0.175,rady=0.15,lwd=2,lab="Zooplankton",font=2)
textempty(mid=c(0.85,0.4),lab="Respiration",font=3)
textempty(mid=c(0.15,0.4),lab="Predation",font=3)
textempty(mid=c(0.5,0.8),lab="Ingestion",font=3)
textempty(mid=c(0.5,0.25),lab="Defaecation",font=3)
box(col="grey")
subtitle()

########################################
# fig 1.3: scientific method
########################################

xx  <- 0.2
yy  <- 0.03
openplotmat()
textrect(c(0.20,0.97),xx,yy,lab="Real world",box.col="lightgrey",shadow.size=0)
textrect(c(0.78,0.97),xx,yy,lab="Conceptual world",box.col="lightgrey",shadow.size=0)
abline(v=0.5)
straightarrow (from=c(0.20,0.8) ,to=c(0.55,0.8) ,lwd=2,arr.pos=1)
straightarrow (from=c(0.55,0.15),to=c(0.18,0.15),lwd=2,arr.pos=1)
straightarrow (from=c(0.78,0.7) ,to=c(0.78,0.62),lwd=2,arr.pos=1)
straightarrow (from=c(0.78,0.43),to=c(0.78,0.35),lwd=2,arr.pos=1)

textrect(c(0.20,0.5),xx,yy,lab="Phenomena",shadow.size=0)
textrect(c(0.78,0.8),xx,yy,lab="Observations",shadow.size=0)
textrect(c(0.78,0.5),xx,yy,lab="Models",shadow.size=0)
textrect(c(0.78,0.15),xx,yy*4,lab=c("analysis","interpolation","quantification","prediction"),shadow.size=0)
subtitle()

########################################
# fig 1.4: models as interpolation tools
########################################

par(mar=c(0,0,0,0))
x   <- seq(1,186,by=2)
dat <- c(0,sin(2*pi*x/365)*(1+(runif(length(x))-0.5)*0.5)  +(runif(length(x))-0.5)*0.25,0)
x   <- 2*c(0,x,186)
dat <- pmax(0,dat)
#dat <- c(0,5,4,9,10,15,22,20,25,27,24,30,29,25,18,22,15,20,22,20,23,22,18,20,12,15,10,5,7,4,2,0)
#x   <- 1:length(dat)
plot(x,dat,type="l",frame.plot=TRUE,xlab="-",ylab="")
polygon(c(x,1),c(dat,0),col=grey(0.95))

by <- 5
dat2 <- dat[seq(1,length(dat),by=by)]
x2   <- x[seq(1,length(dat),by=by)]
lines(x2,dat2,pch=16,cex=1.5,type="b",lwd=2)


dat3 <- c(dat[1],dat[length(dat)*2/3],dat[length(dat)])
x3   <- c(x[1],x[length(dat)*2/3],x[length(dat)])
lines(x3,dat3,lty=2,lwd=2,pch=18,type="b",cex=1.5)
subtitle()

########################################
# fig 1.5: black box versus models
########################################
par(mfrow=c(2,2))
par(mar=c(1,1,1,1))
openplotmat()
pos <- coordinates (c(3,4,5),hor=FALSE) 
 pos[,2]<-pos[,2]*0.5+0.4
A <- matrix(nr=12,nc=12,1)
A[row(A)<col(A)]<-0
A[row(A)>col(A)+5]<-0
curve<-A
cc <- 0.4
curve[]<-0
curve[3,1]<--cc
curve[6,4]<- cc
curve[7,5]<- cc
curve[7,4]<- cc
curve[10,8]<- curve[11,8]<- curve[12,8]<- cc
curve[11,9]<- curve[12,9]<- curve[12,10]<-cc
A[1,2]<-A[2,3]<-A[4,5]<-A[5,6]<-A[6,7]<-1
A[11,12]<-A[10,11]<-A[9,10]<-A[8,9]<-A[7,8]<-1
A[8,3]     <-0
textrect(c(0.5,0.6),0.5,0.3)
plotmat(A,pos=pos,box.type="rect",box.size=0.08,curve=curve,arr.pos=0.7,
self.shiftx=-0.1,cex=0,name="",box.prop=0.4,arr.len=0.25,self.cex=0.5,relsize=01,add=TRUE)

textrect(c(0.75,0.35),0.25,0.05,box.col="black",lab="black box",col="white",cex=2)

bigarr <- function(x0,y0,x1,y1,col1="darkgrey",col2="black",dx=0.005,dy=0.01)
{

Arrows(x0-dx,y0+dy,x1-dx,y1+dy,arr.col=col2,lwd=25,lcol=col2,
       arr.len=0.7,arr.width=1,arr.type="triangle",arr.adj=0.)
Arrows(x0,y0,x0,y1,arr.col=col1,lwd=25,lcol=col1,
       arr.len=0.7,arr.width=1,arr.type="triangle",arr.adj=0)
}
#bigarr(0.5,0.97,0.5,0.9)
#bigarr(0.5,0.28,0.5,0.21)
textplain(c(0.5,0.225),lab=expression(r^2>0.9),cex=2.5)
textplain(c(0.5,0.1),lab="understanding ??",cex=2.5)

# conceptual model
openplotmat()
textrect(c(0.5,0.6),0.5,0.3)

names <- c("PHYTO","DIN","ZOO","DETRITUS")
M <- matrix(nrow=4,ncol=4,byrow=TRUE,data=c(
#   p n z  d
    0,1,0, 0, #p
    0,0,4, 6,#n
    2,0,0, 0, #z
    0,0,5,0  #d
    ))

pos <- coordinates(c(1,2,1))
pos[,2]<-pos[,2]*0.6+.3

pp<-plotmat(M,pos=pos,curve=0,name=names,lwd=1,my=0.0,cex.txt=0.8,prefix="f",
            box.lwd=2,box.size=0.08,box.type="square",box.prop=0.5,box.cex=0.8,
            arr.type="triangle",arr.pos=0.6,shadow.size=0.01,add=TRUE)

# extra arrows: flow 5 to Detritus and flow 2 to detritus
phyto   <-pp$comp[names=="PHYTO"]
zoo     <-pp$comp[names=="ZOO"]
nh3     <-pp$comp[names=="DIN"]
detritus<-pp$comp[names=="DETRITUS"]

# flow2->detritus
m2 <- 0.5*(zoo+phyto)
m1 <- detritus
m1[1] <-m1[1]+ pp$radii[3,1]*0.2
m1[2] <-m1[2] + pp$radii[3,2]
mid<-straightarrow (to=m1,from=m2,arr.type="triangle",arr.pos=0.7,lwd=1)
text(mid[1]-0.01,mid[2]+0.03,"f3",cex=0.8)

# solar radiation
m1 <- 0.5*(nh3+phyto)
m2 <- c(0.25,0.8)
segments (m1[1],m1[2],m2[1],m2[2],lwd=1,lty=2)
text(m2[1]-0.01,m2[2]+0.03,"solar radiation",adj=c(0.5,0.5))

# chlorophyll
m1 <- phyto
m1[1] <-m1[1]+ pp$radii[1,1]
m2 <- m1
m2[1]<-m2[1]+0.25
segments (m1[1],m1[2],m2[1],m2[2],lwd=1)
textellipse(m2,pp$radii[1,1],pp$radii[1,2],lwd=1,shadow.size=0,
 lab="Chlorophyll",cex=0.7)
textrect(c(0.2,0.35),0.2,0.05,lab="model",cex=2)
# DOES NOT FIT IN THIS SIZE WINDOW
#bigarr(0.5,0.97,0.5,0.9)
#bigarr(0.5,0.28,0.5,0.21)
textplain(c(0.5,0.225),lab=expression(r^2<0.9),cex=2.5)
textplain(c(0.5,0.1),lab="!",cex=4.5)
subtitle()

###################################
# Fig. 1.6. Oxygen model fitting
###################################
par(mar=c(5.1,4.1,4.1,2.1))
par(mfrow=c(1,1))
BW   <- 200          # mmol/m3   bottom water conc
Ds   <- 1            # cm2/d     diffusion coeff


O2conc    <- function(Ds,k,BW,x)   BW*exp(-sqrt(k/Ds)*x)
O2flux    <- function(Ds,k,BW)     Ds*sqrt(k/Ds)*BW


x   <- seq(0,4,length=200)
x2   <- seq(0,1,length=200)

# plot flux and penetration depth to screen:
O2flux(Ds,c(1,6,100),BW)/100000*365   # mol/m2/yr

cflux <- function(k,target)Ds*sqrt(k/Ds)*BW/100000*365*12-target

k1  <- uniroot(cflux,c(0,100),target=1  )$root   # gC/m2/yr
k10 <- uniroot(cflux,c(0,100),target=10 )$root   # gC/m2/yr
k100<- uniroot(cflux,c(0,1000),target=100)$root   # gC/m2/yr

O210 <-O2conc(Ds,k=k10,BW,x)*(1+rnorm(200,0,0.05))
plot (O210,x,ylim=rev(range(x)),xlim=c(0,BW),
      xlab=expression("mumol l"^{~-1}),ylab="sediment depth, cm",main=expression(O[2]),col="grey")
lines(O2conc(Ds,k=k10 ,BW,x),x,lty=1,lwd=1,col="black")
legend("bottomright",pch=c(1,NA),lty=c(NA,1),c("measured concentration", "model output"))
subtitle()
########################################
# Fig. 1.7. MODELLING STEPS and INGREDIENTS
########################################

par(mar=c(0,0,0,0))
openplotmat()
elpos<-coordinates (c(1,1,1,1,1,1,1,1),mx=-0.1)
segmentarrow(elpos[7,],elpos[2,],arr.pos=0.15,dd=0.3,arr.side=3)
segmentarrow(elpos[7,],elpos[3,],arr.pos=0.15,dd=0.3,arr.side=3)
segmentarrow(elpos[7,],elpos[4,],arr.pos=0.15,dd=0.3,arr.side=3)

pin   <- par ("pin")        # size of plotting region, inches
xx  <- 0.2
yy  <- xx*pin[1]/pin[2]*0.15  # used to make circles round

sx    <- rep(xx,8)
sx[7] <- 0.05

sy    <- rep(yy,8)
sy[6] <-yy*1.5
sy[7] <- sx[7]*pin[1]/pin[2]

for (i in c(1:7)) straightarrow (from=elpos[i,],to=elpos[i+1,],lwd=2,arr.pos=0.5)
lab <- c("Problem","Conceptual model","Mathematical model","Parameterisation",
         "Mathematical solution","","OK?","Prediction, Analysis")

for (i in c(1:6,8)) textround(elpos[i,],sx[i],sy[i],lab=lab[i])

textround(elpos[6,],xx,yy*2,lab=c("Calibration,sensitivity","Verification,validation"))
textdiamond(elpos[7,],sx[7],sy[7],lab=lab[7])

textplain(c(0.7,elpos[2,2]),yy*2,lab=c("main components","relationships"),font=3,adj=c(0,0.5))
textplain(c(0.7,elpos[3,2]),yy ,"general theory",adj=c(0,0.5),font=3)
textplain(c(0.7,elpos[4,2]),yy*2,lab=c("literature","measurements"),font=3,adj=c(0,0.5))
textplain(c(0.7,elpos[6,2]),yy*2,lab=c("field data","lab measurements"),font=3,adj=c(0,0.5))

subtitle()

par(ask=opar$ask)
par(mar=opar$mar)
par(oma=opar$oma)
par(mfrow=opar$mfrow)

