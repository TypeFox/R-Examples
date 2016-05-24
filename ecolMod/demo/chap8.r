##############################
## Soetaert and Herman      ##
## ecological modelling     ##
## Figures from chapter 8   ##
## Equilibrium analysis     ##
##############################

opar <- par()
par(ask=TRUE)
figNr <- 1
subtitle <- function()
{
 mtext(side=1,outer=TRUE,"Soetaert and Herman - chapter 8  ",cex=0.7,adj=1,line=-1.5)
 mtext(side=1,outer=TRUE,paste("  Fig. 8.",figNr,sep=""),cex=0.7,adj=0,line=-1.5)
 figNr <<- figNr +1
}


###############################################################################
####======================================================================#####
####                                Theory                                #####
####======================================================================#####
###############################################################################

################################################
# Fig. 8.1 Fraction ammonia/total ammonia, vs pH  
################################################

par(mfrow=c(2,2))
par(mar=c(1,1,1,1))
openplotmat()
elpos<-coordinates (c(1,2),hor=FALSE)
tt <- treearrow(from=elpos[1,],to=elpos[2:3,],lty=1,path="V")
text (tt[1,1],tt[1,2]+0.05,expression(k^"+"))
text (tt[2,1],tt[2,2]+0.05,expression(k^"+"))
tt <- treearrow(from=elpos[2:3,],to=elpos[1,],lty=1,path="V")
text (tt[1,1],tt[1,2]+0.05,expression(k^"-"))
names<- c(expression(NH[4]^"+"),expression(NH[3]),expression(H^"+"))
for ( i in 1:3) textrect (elpos[i,],0.08,0.08,lab=names[i],cex=1.5)

#box(col="grey")
par(mar=c(5.1,4.1,4.1,2.1))
writelabel("A",line=1,at=0.1)

par(mar=c(5.1,4.1,4.1,2.1))
require(seacarb)
KN           <- Kn(0,30,0)         # mol/kg , equilibrium ct at 30dgC
curve (KN/(KN+10^(-x)),from=6,to=12,lwd=2,xlab="pH",ylab="-",main="fraction ammonia")
KN           <- Kn(0,0,0)          # mol/kg , equilibrium ct at 0 dgC
curve (KN/(KN+10^(-x)),lwd=1,lty=2,add=TRUE)
legend("bottomright",c("30 dgC","0dgC"), lty=c(1,2),lwd=c(2,1))
writelabel("B")
subtitle()

################################################
# Fig. 8.2 Ammonium equilibrium + slow reaction
################################################

par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
openplotmat()
rect(0.2,0.05,0.8,0.65, angle=45,density=15,col="darkgrey",border=NA)
elpos<-coordinates (c(1,2),hor=FALSE,my=-0.15,mx=0.0,relsize=0.8)
tt <- treearrow(from=elpos[1,],to=elpos[2:3,],lty=1,path="V")
text (tt[1,1],tt[1,2]+0.05,expression(k^"+"))
text (tt[2,1],tt[2,2]+0.05,expression(k^"+"))
tt <- treearrow(from=elpos[2:3,],to=elpos[1,],lty=1,path="V")
text (tt[1,1],tt[1,2]+0.05,expression(k^"-"))
tt<-bentarrow(from=elpos[2,],to=elpos[2,]+c(0.2,-0.1),arr.pos=1)
text(tt[1,1]+0.025,tt[1,2],expression(lambda),cex=1.2)
names<- c(expression(NH[4]^"+"),expression(NH[3]),expression(H^"+"))
for ( i in 1:3) textrect (elpos[i,],0.08,0.08,lab=names[i],cex=1.5)
box(col="grey")

par(new=TRUE)
par(fig=c(0,0.35,0.65,1.0))
par(mar=c(1,1,1,1))
openplotmat()
elpos<-c(0.5,0.5)
tt<-bentarrow(from=elpos,to=elpos+c(0.4,-0.2),arr.pos=1)
text(tt[1,1]+0.05,tt[1,2],expression(lambda),cex=1.2)
textrect (elpos,0.2,0.2,lab=expression(sum(NH[x])),cex=1.5)
box(col="grey")
par(new=FALSE)
par(fig=c(0,1,0.0,1))

subtitle()


################################################
# Fig. 8.3. Equilibrium Monod
################################################
par(mar=c(1,1,1,1))
openplotmat()
rect(0.075,0.05,0.575,0.45, angle=45,density=15,col="darkgrey",border=NA)

elpos<-coordinates (c(4,2,4),hor=FALSE)
elpos<-elpos[-c(3,4,7,10),]
treearrow(from=elpos[1:2,],to=elpos[3,],lty=1,path="V")
treearrow(from=elpos[3,],to=elpos[1:2,],lty=1,path="V")
treearrow(from=elpos[3:4,],to=elpos[5:6,],lty=1,path="V")
names<- c("A","E","EA","B","S","E")
for ( i in 1:6) textrect (elpos[i,],0.06,0.06,lab=names[i],cex=1.5)
text(0.4,0.28,expression(k^{"+"}))
text(0.3,0.15,expression(k^{"-"}))
text(0.3,0.4,expression(k^{"-"}))
text(0.735,0.4,"r")
text(0.735,0.65,"r")
box(col="grey")

par(new=TRUE)
par(fig=c(0,0.4,0.6,1.0))
par(mar=c(1,1,1,1))
openplotmat()
elpos<-coordinates (c(2,1),hor=FALSE)
treearrow(from=elpos[1:2,],to=elpos[3,],lty=1,path="V")
names<- c("A","B","S")
for ( i in 1:3) textrect (elpos[i,],0.09,0.09,lab=names[i],cex=1.5)
text(0.55,0.55,expression(r[f]))
box(col="grey")
par(fig=c(0,1,0.0,1))

subtitle()

################################################
# Fig. 8.4. Transport in porous media
################################################

par(mar=c(0,0,0,0))
emptyplot(c(0,1),c(0,1),asp=TRUE)
rect(0.1,0.1,0.9,0.9,col=grey(0.95))

wx   <- 0.14
wy   <- 0.14
filledmultigonal(mid=c(0.5,0.5),rx=wx,ry=wy,nr=6,col=grey(0.4))

text(0.5,0.5,"S",font=2,cex=2)
text(0.3,0.7,"C",font=2,cex=1.5)
text(0.75,0.55,"C",font=2,cex=1.5)
text(0.5,0.25,"C",font=2,cex=1.5)

Arrows(0.45,0.55,0.35,0.65,code=3)
Arrows(0.55,0.51,0.7,0.535,code=3)
Arrows(0.5,0.45,0.5,0.3,code=3)
bentarrow(c(0.775,0.55),c(0.85,0.5),arr.adj=0,arr.len=0.25,lwd=1)
bentarrow(c(0.275,0.7),c(0.2,0.65),arr.adj=0,arr.len=0.25,lwd=1)
bentarrow(c(0.525,0.25),c(0.6,0.2),arr.adj=0,arr.len=0.25,lwd=1)

plotellipse(mid=c(0.32,0.85),rx=0.05,ry=0.05,arrow=TRUE,arr.pos=1,angle=55,from=0,to=2*pi,arr.length=0.3,lwd=1)
plotellipse(mid=c(0.2,0.45),rx=0.05,ry=0.05,arrow=TRUE,arr.pos=1,angle=95,from=0,to=2*pi,arr.length=0.3,lwd=1)
plotellipse(mid=c(0.75,0.75),rx=0.05,ry=0.05,arrow=TRUE,arr.pos=1,angle=255,from=0,to=2*pi,arr.length=0.3,lwd=1)
plotellipse(mid=c(0.3,0.25),rx=0.05,ry=0.05,arrow=TRUE,arr.pos=1,angle=95,from=0,to=2*pi,arr.length=0.3,lwd=1)
plotellipse(mid=c(0.82,0.35),rx=0.05,ry=0.05,arrow=TRUE,arr.pos=1,angle=255,from=0,to=2*pi,arr.length=0.3,lwd=1)
subtitle()


###############################################################################
####======================================================================#####
####                             R case study                             #####
####======================================================================#####
###############################################################################


##################################
# 8.5 pH in algal culture
##################################
# the dissociation constants
Salinity     <- 0
Temperature  <- 20
WDepth       <- 0

k1           <- K1(Salinity,Temperature,WDepth)    # Carbonate k1
k2           <- K2(Salinity,Temperature,WDepth)    # Carbonate k2

par (mar=c(5.1,4.1,4.1,2.1))

pHfunction <- function(pH, k1,k2, DIC, Alkalinity )
{
   H    <- 10^(-pH)
   HCO3 <- H*k1  /(H*(k1+H) + k1*k2)*DIC
   CO3  <- k1*k2 /(H*(k1+H) + k1*k2)*DIC

   EstimatedAlk  <- (- H) *1.e6  + HCO3 + 2*CO3

   return(EstimatedAlk  - Alkalinity)

}

Alkalinity <- 2200
DIC        <- 2100

sol <- uniroot(pHfunction,lower=0,upper=12, tol=1.e-20,
                k1=k1, k2=k2, DIC=DIC, Alkalinity=Alkalinity)
sol$root


#========================================================================
# pH changes due to primary production
#========================================================================
# load package with the integration routine:

require(deSolve)

#----------------------#
# the model equations: #
#----------------------#
model<-function(time,state,parameters)
  {
with(as.list(c(state,parameters)),{

    PAR    <- 0.
    if(time%%24 < dayLength) PAR <- parDay

    Growth <- maxGrowth*DIN/(DIN+ksDIN)*PAR/(PAR+ksPAR)*ALGAE -
                respRate * ALGAE

    dDIN   <- -Growth                   # DIN is consumed
    dDIC   <- -Growth * CNratio         # DIC is consumed ~ CN ratio
    dALGAE <- Growth                    # algae increase by growth
    dALKALINITY <- Growth               #alkalinity production if nitrate


    # estimate the pH
    pH  <- uniroot(pHfunction,lower=0,upper=12,tol=1.e-20,
           k1=k1,k2=k2, DIC=DIC,Alkalinity=ALKALINITY)$root

   list(c(dDIN,dALGAE,dALKALINITY,dDIC),c(PAR=PAR,pH=pH)  )

    })
 }

#-----------------------#
# the model parameters: #
#-----------------------#

parameters<-c(maxGrowth   =0.125,      #molN/molN/hr  Maximal growth rate
              ksPAR       =100,        #muEinst/m2/s  Half-saturation ct for light-limited growth
              ksDIN       =1.0,        #mmolN/m3      Half-saturation ct of N uptake Phytoplankton
              respRate    =0.001,      #/h            Respiration rate
              CNratio     =6.5,        #molC/molN     carbon:Nitrogen ratio
              parDay      =250.,       #muEinst/m2/s  PAR during the light phase
              dayLength   =12.         #hours         Length of illuminated period (in one day)
              )
Salinity     <- 0
Temperature  <- 20
WDepth       <- 0

# the dissociation constants
k1           <- K1(Salinity,Temperature,WDepth)    # Carbonate k1
k2           <- K2(Salinity,Temperature,WDepth)    # Carbonate k2

#-------------------------#
# the initial conditions: #
#-------------------------#

state     <-c(DIN        =30,     #mmolN/m3
              ALGAE      =0.1,    #mmolN/m3
              ALKALINITY =2200,   #mmol/m3
              DIC        =2100)   #mmolC/m3

#----------------------#
# RUNNING the model:   #
#----------------------#

times <-seq(0,24*10,1)

require(deSolve)
out   <-as.data.frame(ode(state,times,model,parameters))


#------------------------#
# PLOTTING model output: #
#------------------------#

par(mfrow=c(2,2), oma=c(0,0,3,0))         # set number of plots (mfrow) and margin size (oma)
plot (out$time,out$ALGAE,type="l",main="Algae",xlab="time, hours",ylab="mmol/m3")
polygon(out$time,out$PAR-10,col="lightgrey",border=NA)
box()
lines (out$time,out$ALGAE  ,lwd=2 )

writelabel("A")
plot (out$time,out$DIN ,type="l",main="DIN"  ,xlab="time, hours",ylab="mmolN/m3", lwd=2)
writelabel("B")
plot (out$time,out$DIC ,type="l",main="DIC"  ,xlab="time, hours",ylab="mmolC/m3", lwd=2)
writelabel("C")
plot (out$time,out$pH,type="l",main="pH"  ,xlab="time, hours",ylab="-", lwd=2)
writelabel("D")
mtext(outer=TRUE,side=3,"Algal growth and pH",cex=1.5)

subtitle()
