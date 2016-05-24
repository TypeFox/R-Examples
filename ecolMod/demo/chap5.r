##############################
## Soetaert and Herman      ##
## ecological modelling     ##
## Figures from chapter 5   ##
## Analytic solutions       ##
##############################
opar <- par()
par(ask=TRUE)
par(mfrow=c(1,1))
figNr <- 1
subtitle <- function()
{
 mtext(side=1,outer=TRUE,"Soetaert and Herman - chapter 5  ",cex=0.7,adj=1,line=-1.5)
 mtext(side=1,outer=TRUE,paste("  Fig. 5.",figNr,sep=""),cex=0.7,adj=0,line=-1.5)
 figNr <<- figNr +1
}

###############################################################################
####======================================================================#####
####                                Theory                                #####
####======================================================================#####
###############################################################################

par(mfrow=c(1,2))
par(mar=c(1,1,1,1),cex=0.7)
openplotmat()
pos <-coordinates(c(2,2,2,1,1), hor=TRUE)
pos <- pos[-c(1,3),]
pos[1,2]<-pos[1,2]+0.04
pos[2,2]<-pos[2,2]+0.025
pos[5,2]<-pos[5,2]-0.05
pos[6,2]<-pos[6,2]-0.03

straightarrow(pos[1,],pos[2,],arr.pos=0.7)
straightarrow(pos[2,],pos[4,],arr.pos=0.7)
treearrow(from=pos[3:4,],to=pos[5,],arr.pos=0.4)
straightarrow(pos[5,],pos[6,],arr.pos=0.6)

textrect(pos[2,],radx=0.1,rady=0.05,lab="Integration",shadow.size=0)
textround(pos[1,],0.1,0.05,lab="Differential equation",shadow.size=0.01)
textround(pos[3,],0.1,0.05,lab= c("Initital and/or","boundary conditions"),shadow.size=0.01)
textround(pos[4,],0.1,0.05,lab="General solution",shadow.size=0.01)
textrect(pos[5,],radx=0.15,rady=0.06,lab=c("Solving integration","constants"),shadow.size=0)

textround(pos[6,],0.1,0.05,lab="Particular solution",shadow.size=0.01)
box(col="grey")
writelabel("A",line=-3,at=0.1)

par(mar=c(1,1,1,1),cex=0.8)
pos <-coordinates(c(2,2,2,1,1), hor=TRUE)
pos <- pos[-c(1,3),]
pos[1,2]<-pos[1,2]+0.04
pos[2,2]<-pos[2,2]+0.025
pos[5,2]<-pos[5,2]-0.05
pos[6,2]<-pos[6,2]-0.03

openplotmat()
straightarrow(pos[1,],pos[2,],arr.pos=0.7)
straightarrow(pos[2,],pos[4,],arr.pos=0.7)
treearrow(from=pos[3:4,],to=pos[5,],arr.pos=0.4)
straightarrow(pos[5,],pos[6,],arr.pos=0.6)

textrect(pos[2,],radx=0.1,rady=0.05,lab="Integration",shadow.size=0)
textround(pos[1,],0.1,0.05,lab=expression(frac(dC,dt)==kC),shadow.size=0.01)
textround(pos[3,],0.1,0.05,lab=expression(C["t=0"]==C[0]),shadow.size=0.01)
textround(pos[4,],0.1,0.05,lab=expression(C==A*e^{kt}),shadow.size=0.01)
textrect(pos[5,],radx=0.15,rady=0.06,lab=c(expression(C[0]==A*e^{k*0}),expression(C[0]==A)),shadow.size=0)

textround(pos[6,],0.1,0.05,lab=expression(C==C[0]*e^{kt}),shadow.size=0.01)
box(col="grey")
writelabel("B",line=-3,at=0.1)
subtitle()

########################################
# Figure 5.2. Shape of exponential function
########################################

par(mfrow=c(2,2),mar=c(4,4,3,1))
curve(exp(0.1*x),0,150,main=expression(y==exp^(0.1*x)),xlab="x",ylab="y",lwd=2)
writelabel("A")
curve(exp(-0.1*x),0,150,main=expression(y==exp^(-0.1*x)),xlab="x",ylab="y",lwd=2)
writelabel("B")
subtitle()

########################################
# Figure 5.3. O2-BOD curve
########################################

#pp <- par(mar=c(0,0,0,0))
par(mfrow=c(1,2))
openplotmat()
rect(0.0,0.0,1,1)
rect(0.01,0.01,0.99,0.7,border=NA,col="lightgrey")
elpos <- matrix(nrow=2,data=c(0.25,0.75,0.4,0.4))
straightarrow(from=elpos[1,]+c(0,0.1),to=elpos[1,]+c(0,0.3),arr.pos=1)
straightarrow(to=elpos[1,]+c(0,0.1),from=elpos[1,]+c(0,0.3),arr.pos=1,arr.adj=1)
text(elpos[1,1],elpos[1,2]+0.4,expression(O[2]^"*"))
bentarrow (from=elpos[1,],to=elpos[1,]+c(0.2,-0.2),arr.pos=1.0,lwd=2)
bentarrow (from=elpos[2,],to=elpos[2,]+c(-0.2,-0.2),arr.pos=1.0,lwd=2)

textrect(elpos[1,],0.1,0.1,lab=expression(O[2]),shadow.size=0.02)
textrect(elpos[2,],0.1,0.1,lab="BOD",shadow.size=0.02)
writelabel("A")

#par(mar=pp)

k       = 0.1             # /day     - reaeration
O2sat   = 300             # mmol/m3  - saturated oxygen concentration
r       = 0.05            # /day     - BOD decay rate
O2_0    = 250             # mmol/m3  - Initial oxygen concentration
BOD_0   = 500             # mmol/m3  - Initial BOD concentration
par(mar=c(5.1,4.1,4.1,2.1))
curve(BOD_0*exp(-r*x),0,100,xlab="time",ylab="mmol O2/m3",lwd=2)
curve(BOD_0*r*(exp(-k*x)-exp(-r*x)) /(k-r)+O2_0*exp(-k*x)+O2sat*(1-exp(-k*x)),0,100,lty=2,lwd=2,add=TRUE)
abline(h=O2sat,lty=3)
legend("topright",c("BOD","O2",expression(O[2]^"*")),lwd=c(2,2,1),lty=c(1,2,3),cex=0.8)

writelabel("B")

subtitle()


###############################################################################
####======================================================================#####
####                            R case studies                            #####
####======================================================================#####
###############################################################################

########################################
# Figure 5.4. 
# The diffusion-reaction equation in 1-D 
# cartesian coordinates  
########################################

Ds  <- 1   # diffusion coefficient
ini <- 1   # initial condition
k   <- 0.05 # growth rate

grow1D <- outer(xx<-seq(-5,5,length=50),tt<-seq(0.1,5,by=0.025),
                FUN = function (x,tt) ini/(2*sqrt(pi*Ds*tt))*exp(k*tt-x^2/(4*Ds*tt)))
par(mfrow=c(1,1))
par(mar=c(5.1,4.1,4.1,2.1))
persp(xx,tt,z=grow1D,theta=150,box=TRUE,axes=TRUE,col=drapecol(grow1D,femmecol(100)),
xlab="space",ylab="time",zlab="Density",border=NA,main="1-D dispersion-reaction")  

subtitle()

########################################
# Figure 5.5. 
# The diffusion-reaction equation in 1-D 
# cylindrical coordinates  
########################################
                         
growplane <- outer(rr<-seq(-5,5,length=50),tt<-seq(0.1,5,by=0.025),
              FUN = function (rr,tt) ini/(4*pi*Ds*tt)*exp(k*tt-rr^2/(4*Ds*tt)))

plotplane <- function(time,rmax=5,...) 
 { 
  val <-outer(xx<-seq(-rmax,rmax,length=50),yy<-xx,
              FUN = function (x,y) 
             {
              r2<-x*x+y*y;
              ini/(4*pi*Ds*time)*exp(k*time-r2/(4*Ds*time))
              } 
             ) 
  persp(xx,yy,z=val,theta=150,box=TRUE,axes=TRUE,col=drapecol(val,femmecol(100)),zlab="Density",border=NA,...)         
 }

par(mfrow=c(2,2),mar=c(3,3,3,3))
plotplane(0.1, main= "0.1 day")
plotplane(1  , main= "1 day")
plotplane(2  , main= "2 days") 
plotplane(5  , main= "5 days")  

subtitle()



##########################################
# Fig. 5.6. 1-dimensional oxygen model in organism
##########################################

# parameters

BW     <- 2            # mmol/m3,       oxygen concentration in surrounding water
Da     <- 0.5          # cm2/d          effective diffusion coefficient in organism
R      <- 0.005        # cm             radius of organism 
Q      <- 250000       # nM/cm3/d       oxygen consumption rate per volume per day

rr        <- seq(-R,R,length=400)

sandwich <- function(Da,Q,BW,R,r)  BW+Q/(2*Da)*(r^2-R^2)
cylinder <- function(Da,Q,BW,R,r)  BW+Q/(4*Da)*(r^2-R^2)
sphere   <- function(Da,Q,BW,R,r)  BW+Q/(6*Da)*(r^2-R^2)
par(mfrow=c(2,2))
 par(mar=c(5.1,4.1,4.1,2.1))
plot  (rr,sandwich(Da,Q,BW,R,rr),lwd=2,lty=1,
       ylim=c(0,BW), type="l",main="oxygen in organism" ,
       xlab="radius,cm",ylab="mumol/l")
lines (rr,cylinder(Da,Q,BW,R,rr),lwd=2,lty=2)
lines (rr,sphere  (Da,Q,BW,R,rr),lwd=2,lty=3)
legend("top",c("sheet","cylinder","sphere"),lty=1:3,lwd=2)
writelabel("A")

plot  (rr,cylinder(Da,Q,BW,R,rr),ylim=c(0,BW), type="l",lwd=2,
       xlab="radius,cm",ylab="O2, mumol/l",main="cylinder" )
lines (rr*0.5,cylinder(Da,Q,BW,R*0.5,rr*0.5),lwd=2,lty=2)
lines (rr*0.25,cylinder(Da,Q,BW,R*0.25,rr*0.25),lwd=2,lty=3)
legend("bottom",legend=c("100mum","50mum","25mum"),title="maximal thickness",lty=1:3,lwd=2)
writelabel("B")

critsand <- function(Da,Q,BW)  sqrt(BW*2*Da/Q)
critcyl  <- function(Da,Q,BW)  sqrt(BW*4*Da/Q)
critsph  <- function(Da,Q,BW)  sqrt(BW*6*Da/Q)

BWseq <- seq(0,50,length=100)
plot(BWseq,10000*sqrt(BWseq*6*Da/Q),type="l",lty=3,lwd=2,
    main="critical thickness", xlab="surrounding oxygen, mumol/l",ylab="mum")
lines(BWseq,10000*sqrt(BWseq*4*Da/Q),lty=2,lwd=2)
lines(BWseq,10000*sqrt(BWseq*2*Da/Q),lty=1,lwd=2)
legend("bottomright",c("sheet","cylinder","sphere"),lty=1:3,lwd=2)

writelabel("C")

emptyplot(c(-1.5,1.5),main="shapes" )
ypos <- c(0.9,1.1)
plotellipse(mid=c(0,ypos[1]),rx=1,ry=0.4,from=-pi,to=0,lwd=1)

plotellipse(mid=c(0,ypos[2]),rx=1,ry=0.4,col="lightgrey",lwd=1)
plotellipse(mid=c(0,ypos[1]),rx=1,ry=0.4,from=0,to=pi,lwd=1,lcol="darkgrey")
segments(-1,ypos[1],-1,ypos[2],lwd=2)
segments( 1,ypos[1], 1,ypos[2],lwd=2)

rseq1 <- seq(0.3,1,length.out=50)
col  <- grey(c(rseq1,rev(rseq1)))
filledcylinder(rx=0.1,ry=0.2,len=2.5,angle=0,col=col,mid=c(0,0),lcol="black",lwd=1)

col2  <- grey(seq(1,0.5,length.out=50))
filledcircle(r1=0.5,col=col2,mid=c(0,-1.),lcol="black",lwd=1)
box()
writelabel ("D")
subtitle()

########################################
# Figure 5.7. Nonlocal exchage model  

nonlocal <- function( 
        k=0.03108,v=0.001,Db=0.1,depo=0.1,L=5,injectflux=0.3,sed
                     )
{
 # Power and exponents
 a  <- (v/Db - sqrt ( (v/Db)^2 + 4*k / Db))/2.
 b  <- (v/Db + sqrt ( (v/Db)^2 + 4*k / Db))/2.
 expaL <- exp(a*L)
 expbL <- exp(b*L)

 A <- matrix(nrow=3,ncol=3,byrow=TRUE,  data=   c(

   v-Db*a         ,v-Db*b         ,0             , 
   expaL          ,expbL          ,-expaL        ,   
(v-Db*a)*expaL  ,(v-Db*b)*expbL ,(Db*a-v)*expaL )
           )
# right hand side

 B  <- c(depo     ,  0            ,- injectflux)

 X  <- solve(A,B)

 s1 <- which (sed<L)
 s2 <- which (sed>=L)

 conc <- vector(length=length(sed))
 conc[s1] <- X[1]* exp(a*sed[s1])+X[2]*exp(b*sed[s1])
 conc[s2] <- X[3]* exp(a*sed[s2])
 return(conc)

}   # end nonlocal


#-------------------------------------------------------------------*
# Model application 1: Pb210
#-------------------------------------------------------------------*

depth <- seq(0,15,by=0.01)              # cm          depth of sediment slices

# First set of parameters
k            <- 0.03108                 # /yr         first-order decay rate
Db           <- 0.1                     # cm2/yr      bioturbation coefficient
w            <- 0.001                   # cm/yr       advection rate
depo         <- 0.1                     # dpm/cm2/yr  Pb210 deposition flux
injecflux    <- 0.3                     # dpm/cm2/yr  Pb210 injection flux
injecdepth   <- 5                       # cm          Pb210 injection depth


Pb <-nonlocal(k,w,Db,depo,injecdepth,injecflux,sed=depth)
par(mfrow=c(1,1))
par(mar=c(5.1,4.1,4.1,2.1))
plot(Pb,depth,ylim=c(15,0),xlab="Pb, dpm/cm3",ylab="cm",type="l",lwd=2,main="75% injected")


subtitle()
par(ask=opar$ask)
par(mar=opar$mar)
par(oma=opar$oma)
par(mfrow=opar$mfrow)

