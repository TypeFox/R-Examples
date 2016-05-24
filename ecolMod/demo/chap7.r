##############################
## Soetaert and Herman      ##
## ecological modelling     ##
## Figures from chapter 7   ##
## Stability analysis       ##
##############################

opar <- par()
par(mfrow=c(1,1))

par(ask=TRUE)
figNr <- 1
subtitle <- function()
{
 mtext(side=1,outer=TRUE,"Soetaert and Herman - chapter 7  ",cex=0.7,adj=1,line=-1.5)
 mtext(side=1,outer=TRUE,paste("  Fig. 7.",figNr,sep=""),cex=0.7,adj=0,line=-1.5)
 figNr <<- figNr +1
}

###############################################################################
####======================================================================#####
####                                Theory                                #####
####======================================================================#####
###############################################################################

##############################################
# Fig. 7.1 local stable/unstable equilibria
##############################################
fun <- function(x) x*sin(x)
par(mfrow=c(1,1))
par (mar=c(0,0,0,0))

curve  (fun(x),-2,7,axes=FALSE,xlab="",ylab="",lwd=2,n=500)

minima <- rbind(c(0,0),c(4.91,fun(4.91)))
rad <- seq(0,0.5,length.out=100)
col1<- greycol(100)
col2<- greycol(100,c(0,0.4))   # less dark 

for (i in 1:2) filledellipse(rx1=0.5,mid=minima[i,]+c(0,0.5),col=col1)

filledellipse(rx1=0.5,mid=c(2.0,fun(2.))+c(0,0.5),col=col2)
len <- 0.3
Arrows(-1.5 ,fun(-1.5)+0.5,-1.0 ,fun(-1.0 )+0.5,arr.length=len)
Arrows(2.5 ,fun(2.5 )+0.5,3 ,fun(3 )+0.5,arr.length=len)
Arrows(1.5 ,fun(1.5 )+0.5,1 ,fun(1 )+0.5,arr.length=len)

Arrows(3.8 ,fun(3.8 )+0.5,4.2 ,fun(4.2 )+0.5,arr.length=len)
Arrows(6. ,fun(6. )+1,5.5,fun(5.5)+1,arr.length=len)

subtitle()

##############################################
# Fig. 7.2 Stable/unstable equilibrium - the model equation
##############################################

# the model : density-dependent growth and Monod-type mortality rate
model <- function(N) r*N*(1-N/K)-Q*N/(N+ks)

# parameter values
K  <- 10           # carrying capacity
Q  <- 0.1          # harvesting rate         
ks <- 1            # half-saturation concentration 
r  <-  0.05        # rate of increase
rcrit <- 4*Q*K/(ks+K)^2       # critical rate of increase below which the population becomes extinct

# plot the function for a certain value of r:
plotfun <- function (rin=0.038,lab="")
{
  r <<- rin
  curve(model(x),0,10,n=500,xlab="N",ylab="dN/dt",lwd=2) 
  abline(h=0,lty=2)                # add 0-axis
  legend("bottomleft",paste("r=",r))
  equi  <- equilibrium()
  DrawArrow(equi,dy=0)
  writelabel(lab)
}   

# function to estimate the equilibrium points
equilibrium <- function( )
  {
  DD    <- (ks+K)^2-4*Q*K/r
  if (DD >=0) 
  {
    SD    <- sqrt(DD) 
    eq    <- unique(c(0,0.5*((K-ks)+SD),0.5*((K-ks)-SD)))
  } else eq <- 0
  eq [eq>=0]
}

# function to draw arrows

DrawArrow <- function(equi,d1=1,d2=0.25, dy = 0.01)
# equi: the equilibrium points
# d1:largest distance of arrow to equilibrium point
# d2:smallest distance
# dy: distance of text from equilibrium points

{
  tiny      <- 1e-8   # perturbation factor

  for (ee in equi)
  {
   # arrows: code ~ position of arrow head at end or begin
   ifelse (sign(model(ee-tiny)) == 1, code1 <-2, code1 <- 1)
   arrows(ee-d1,0,ee-d2,0,code=code1,length=0.1,lwd=2)

   ifelse(sign(model(ee+tiny)) == 1, code2 <- 1, code2 <- 2)

   code <- code1
   if (code != code2) code <- 3   # a saddle point...
   arrows(ee+d1,0,ee+d2,0,code=code2,length=0.1,lwd=2)
   points(ee,0,cex=2,bg=c("white","darkgrey","black")[code],pch=21)
   }
   # name the equilibrium points
  if (dy >0) for (i in 1:length(equi)) text(equi[i],dy,paste("N",i-1,sep=""))
}    ## end DrawArrow

##############################################
# Fig. 7.2 : curve, arrows, equilibrium points
##############################################

par(mar=c(5.1,4.1,4.1,2.1))
r <- 0.038
curve(model(x),0,8,n=500,xlab="N",ylab="dN/dt",lwd=2) 
text(4.5,0.005,"+",cex=2)
text(1,-0.0075,"-",cex=2)
text(8,-0.0075,"-",cex=2)
abline(h=0,lty=2)                # add 0-axis
equi <- equilibrium()
DrawArrow(equi,dy=0.003)
legend( "bottomleft",c("unstable","stable"),pch=21,pt.cex=2,pt.bg=c("white","darkgrey","black"))
subtitle()

#################################################
# Fig. 7.3: 4 different cases
#################################################

par(mfrow=c(2,2),mar=c(4,4,3,2))
plotfun(0.02,"A")
plotfun(rcrit,"B")
plotfun(0.05,"C")
plotfun(0.15,"D")
legend( "topright",c("unstable","stable","saddle"),pch=21,pt.cex=2,pt.bg=c("white","darkgrey","black"))

# formal analysis of equilibriums - pos eigenvalues = unstable, neg=stable
tiny      <- 1e-8   # perturbation factor
eigenvalue <- model(equi+tiny)/tiny
sign(eigenvalue)
subtitle()

#################################################
# Fig. 7.4: multiple stable states - bifurcation
#################################################
par(mfrow=c(1,1))
par(mar=c(5.1,4.1,4.1,6.1))
rSeq <- seq(0,0.15,by=0.0001)
plot(0,xlim=range(rSeq),ylim=c(0,10),type="n",xlab=expression(r[i]),ylab="N*")

tiny <- 1e-8
for (r in rSeq) {
equi <-equilibrium( )
eig <-sign( (model(equi+tiny)/tiny)) # eigenvalue:neg=stable,pos=unstable,0=saddle
points(rep(r,length(equi)),equi,pch=21,
col=c(grey(0.9),"black",grey(0.4))[eig+2],
#bg =c("darkblue","black","lightblue")[eig+2])
bg  =c(grey(0.9),"black",grey(0.4))[eig+2])
}

legend("topleft",pch=21,pt.cex=2,c("stable","unstable"),
pt.bg=c(grey(0.9),grey(0.4)))
abline(v=0.02,lty=2)      ; text(0.02,-0.2,"A",cex=1.2)
abline(v=rcrit,lty=2)     ; text(rcrit,-0.2,"B",cex=1.2)
abline(v=0.05,lty=2)      ; text(0.05,-0.2,"C",cex=1.2)
abline(v=0.15,lty=2)      ; text(0.15,-0.2,"D",cex=1.2)
mtext(side=3,outer=TRUE,"Bifurcation diagram",cex=1.5,line=-2)


# Hysteresis
par(new=TRUE)
par(mar=c(14.1,24.1,13.1,1.1))
#par(mfrow=c(1,2))
# and now less dark
xlim=range(rSeq)
ylim=c(0,10)    ;yr<-range(ylim)
plot(0,xlim=xlim,ylim=ylim,type="n",xlab="",ylab="",axes=FALSE,frame.plot=TRUE)
rect(xlim[1],ylim[1]-yr*0.035,xlim[2],ylim[2]+yr*0.035,col="white",border=NA)
EQ <- NULL
tiny <- 1e-8
for (r in rSeq) {
equi <-equilibrium( )
eig <-sign( (model(equi+tiny)/tiny)) # eigenvalue:neg=stable,pos=unstable,0=saddle
points(rep(r,length(equi)),equi,pch=21,cex=0.5,
col =c(grey(0.92),"black",grey(0.6))[eig+2],
bg  =c(grey(0.92),"black",grey(0.6))[eig+2])
#EQ will contain the path of the arrow...
ij <- which.max(equi)
equi <- equi[ij]
eig  <- eig[ij]
if (min(eig)<0 && equi>0 && r<=0.1) EQ <- rbind(EQ,c(r,equi[eig<0]))

}

#rSeq <- seq(0,0.10,by=0.0001)
#plot(0,xlim=c(0,0.15),ylim=c(0,10),type="n",xlab="r",ylab="N*")
Arrows(0.02,0.15,0.1,0.15,arr.adj=1,lwd=2,arr.type="triangle")
text(0.015,0.15,"a",cex=1.2)
Arrows(0.1,0.15,0.1,max(EQ[,2]),arr.adj=1,lwd=2,arr.type="triangle")
text(0.107,0.5,"b",cex=1.2)
EQ[,2]<-EQ[,2]+0.2
EQ[,1]<-EQ[,1]-0.002
lines(EQ,lwd=2)
text(0.105,9.5,"c",cex=1.2)
text(EQ[1,1]-0.005,EQ[1,2],"d",cex=1.2)
text(0.04,0.7,"e",cex=1.2)
Arrows(EQ[1,1],EQ[1,2],EQ[1,1],0.15,arr.adj=1,lwd=2,arr.type="triangle")
subtitle()
par(new=FALSE)

#################################################
# Fig. 7.5: hysteresis
#################################################

#Dynamic run
dyna <- function(t,N,pars)
{
r  <<- rs-cos(t/100*pi)*rs*0.9
list(c(model(N)),c(rr=r))
}


rs <- 0.15
times <- seq(0,200,0.5)
N     <- 0.01
out <- as.data.frame(ode(N,times,dyna,parms=0))
N <- out[nrow(out),2]
out <- as.data.frame(ode(N,times,dyna,parms=0))

par(mar=c(5.1,4.1,4.1,2.1))
par(mfrow=c(1,1))
plot(out$rr,out[,2],type="l",lwd=2,xlab="r",ylab="N",main="Hysteresis")
p1 <- 100
p2 <- p1+1
arrows(out[p1,3],out[p1,2],out[p2,3],out[p2,2],length=0.1,lwd=2)
p1 <- 350
p2 <- p1+1
arrows(out[p1,3],out[p1,2],out[p2,3],out[p2,2],length=0.1,lwd=2)
r <- rs
subtitle()

#################################################
# Fig. 7.6: Phase-plane graph
#################################################


par(mgp=c(1,1,0),mar=c(3,3,3,1))

plot(3.2,4.5,xlab="state variable 1", ylab="state variable 2", axes=FALSE,
     main="Isoclines",xlim=c(0,10),ylim=c(0,10), cex.lab=1.3)
box()
usr <- par("usr")
rect(usr[1], usr[3], usr[2], usr[4], col=grey(0.95), border="black")
segments(1,10,5,0,lwd=2)
segments(0,7,9,0,lwd=2,lty=2)

points(3.2,4.5,pch=15,cex=2,)

arrows(1,2,2,3,lwd=2,length=0.1)
arrows(5.5,6.5,4.5,5.5,lwd=2,length=0.1)
arrows(1,7.7,2,6.3,lwd=2,length=0.1)
arrows(5.5,1.6,4.5,2.8,lwd=2,length=0.1)

text(4,8.5,"0-isocline of SV1")
text(8,3,"0-isocline of SV2")
text(3.5,4.8,"Equilibrium",adj=0)

subtitle()



############################################
# the Lotka-Volterra predator-prey equations:
############################################


  LVPredPrey <- function(Time,State,Pars)   # Estimates rate of change
    {
    with(as.list(c(State,Pars)),
    {
    Ingestion     <- rIngestion * PREY*PREDATOR
    GrowthPrey    <- rGrowth * PREY*(1-PREY/capacity)
    MortPredator  <- rMortality * PREDATOR

    dPREY         <- GrowthPrey - Ingestion
    dPREDATOR     <- Ingestion*assEff -MortPredator
    return(list(c( dPREY, dPREDATOR)))
    })

    } #END LVPredPrey

#---------------------------------------------#
# model
#---------------------------------------------#
predprey<-function(PREY,
                   PREDATOR,
                   rIngestion =0.2,    # /day, rate of ingestion
                   rGrowth    =1.0,    # /day, growth rate of prey
                   rMortality =0.2 ,   # /day, mortality rate of prey
                   assEff     =0.5,    # -, assimilation efficiency
                   capacity   =10  )   # mmol/m3, carrying capacity

{

  Time    <- seq(0,200,by=0.1)
  pars    <- c(rIngestion =rIngestion, rGrowth =rGrowth,   
               rMortality =rMortality,assEff =assEff,capacity =capacity)
  out     <- as.data.frame(ode(c(PREY=PREY,PREDATOR=PREDATOR),Time,LVPredPrey,pars))  # ode is integrator
  return(out)
}  # END predprey

#---------------------------------------------#
# isoclines figure
#---------------------------------------------#
isopredprey <- function(rIngestion =0.2,    # /day, rate of ingestion
                        rGrowth    =1.0,    # /day, growth rate of prey
                        rMortality =0.2 ,   # /day, mortality rate of prey
                        assEff     =0.5,    # -, assimilation efficiency
                        capacity   =10 ,... )   # mmol/m3, carrying capacity
{
  Prey1 <- c(0,capacity)
  Pred1 <- rGrowth /rIngestion *(1-Prey1/capacity)
  Prey2 <- rep(rMortality /(rIngestion*assEff),times=2)
  Pred2 <- ylim
  ctprey <- Prey2[1]
  ctpred <- rGrowth /rIngestion *(1-ctprey/capacity)
  EquiX  <- c(ctprey)
  EquiY  <- c(ctpred)
  plot  (Prey1,Pred1, type="n",   # plot first isocline
        main="", axes=FALSE,
        xlab="Prey",ylab="Predator",
        xlim=xlim,ylim=ylim,...)

  usr <- par("usr")
  #rect(usr[1], usr[3], usr[2], usr[4], col=grey(0.95), border=NA)

  lines (Prey1,Pred1,lwd=2)            # plot first isocline
  lines (Prey2,Pred2,lwd=2,lty=2)       # second isocline

  if (length(EquiX > 0)) points(EquiX,EquiY,pch=15,cex=2.0)
  points(0,0,pch=15)
  points(capacity,0,pch=15)

 }


###############################################
# Fig. 7.7: types of equilibriums , 4 on page
###############################################
insetplot <- function(yvar)
{

usr <- par("usr")

xl <- usr[1]
xu <- usr[2]*0.99
yl <- usr[3]
yu <- usr[4]*0.99

rr <- range(yvar)
N <- length(yvar)

p1 <- 0.75

xx <- seq(from=(xl+(xu-xl)*p1 ),to=xu,length.out=N+1)[-1]
yy <- yl+(yu-yl)*p1 + (yu-yl)*(1-p1)*(yvar-rr[1])/(diff(rr))
rect(xl+(xu-xl)*p1,yl+(yu-yl)*p1,xu,yu,col=grey(0.9),lty=2)

lines(xx,yy)
}

par(mgp=c(1,1,0),mfrow=c(2,2),mar=c(4,3,4,2))

xlim          <- c(0,10)            # axes limits
ylim          <- c(0,10)
#---------------------------------------------#
# Stable focus
#---------------------------------------------#

isopredprey( )
box()
title("Stable focal point")
# 4 runs with different initial conditions
seqn <- predprey(0.5,1)
lines(seqn[,2:3],type="l")
arrows(seqn[10,2],seqn[10,3],seqn[11,2],seqn[11,3],length=0.1,lwd=2)

seqn1 <- predprey(5,6)
lines(seqn1[,2:3],type="l",lty=2)
arrows(seqn1[10,2],seqn1[10,3],seqn1[11,2],seqn1[11,3],length=0.1,lwd=2)

seqn <- predprey(3,2)
lines(seqn[,2:3],type="l",lty=3)
arrows(seqn[10,2],seqn[10,3],seqn[11,2],seqn[11,3],length=0.1,lwd=2)

seqn <- predprey(1,8)
lines(seqn[,2:3],type="l",lty=1)
arrows(seqn[10,2],seqn[10,3],seqn[11,2],seqn[11,3],length=0.1,lwd=2)

## The inset, with part of the axes
insetplot(seqn1[,2])
writelabel("A")

#---------------------------------------------#
# Predator extinction
#---------------------------------------------#
xlim          <- c(0,4)            # axes limits
ylim          <- c(0,4)

isopredprey( capacity=1)
box()
title("Predator extinction")
# 2 runs with different initial conditions
seqn2 <- predprey(2,2,capacity=1)
lines(seqn2[,2:3],type="l")
arrows(seqn2[10,2],seqn2[10,3],seqn2[11,2],seqn2[11,3],length=0.1,lwd=2)

seqn <- predprey(3,3,capacity=1)
lines(seqn[,2:3],type="l")
arrows(seqn[10,2],seqn[10,3],seqn[11,2],seqn[11,3],length=0.1,lwd=2)

seqn2 <- predprey(2,2,capacity=1)
lines(seqn2[,2:3],type="l")
arrows(seqn2[10,2],seqn2[10,3],seqn2[11,2],seqn2[11,3],length=0.1,lwd=2)

seqn3 <- predprey(0.5,0.5,capacity=1)
lines(seqn3[,2:3],type="l")
arrows(seqn3[10,2],seqn3[10,3],seqn3[11,2],seqn3[11,3],length=0.1,lwd=2)

insetplot(seqn3[1:500,2])
writelabel("B")


#---------------------------------------------#
# Neutral stability
#---------------------------------------------#
xlim          <- c(0,10)            # axes limits
ylim          <- c(0,10)
par(new=FALSE)
isopredprey( capacity=1e6)
box()
title("Neutral stability")
# 2 runs with different initial conditions
seqn2 <- predprey(2,2,capacity=1e6)
lines(seqn2[,2:3],type="l")
arrows(seqn2[10,2],seqn2[10,3],seqn2[11,2],seqn2[11,3],length=0.1,lwd=2)

seqn <- predprey(3,3,capacity=1e6)
lines(seqn[,2:3],type="l")
arrows(seqn[10,2],seqn[10,3],seqn[11,2],seqn[11,3],length=0.1,lwd=2)

insetplot(seqn2[1:1000,2])
writelabel("C")

emptyplot()
legend("center",c("constant prey","constant predator"), lwd=2,lty=c(1,2))

subtitle()


###############################################
# Fig. 7.8: Trajectories to equilibrium
###############################################
par(mgp=opar$mgp)
par(mfrow=c(2,2)) 

############################################
# the Lotka-Volterra competition equations:
############################################

competition<-function(SPEC1      =10,        # initial condition of spec 1
                      SPEC2      =6,
                      rGrowth1   =0.1,       # /day, growth rate of spec1
                      rGrowth2   =0.2,       # /day, growth rate of spec2
                      capacity1  =7 ,        # mmol/m3, carrying capacity spec1
                      capacity2  =5,         # mmol/m3, carrying capacity spec2
                      comp12     =0.5,       # competitive effect of sp 1 on growth of spec 2
                      comp21     =0.4)

{

  model <- function(Time,State,Pars)   # Estimates rate of change
    {
    with(as.list(State),
    {
    dSPEC1  <- rGrowth1  * SPEC1*(1-(SPEC1+comp12*SPEC2)/capacity1)
    dSPEC2  <- rGrowth2  * SPEC2*(1-(SPEC2+comp21*SPEC1)/capacity2)
    return(list(c( dSPEC1, dSPEC2)))
    })
    
    } #END model

  Time    <- seq(0,200,by=0.1)
  out     <- as.data.frame(ode(c(SPEC1=SPEC1,SPEC2=SPEC2),Time,model,0))  # ode is integrator
  return(out)
}

seqn3  <- competition(10,9)                        # phase plane plot
seqn4  <- competition(15,6,capacity1 =17,comp12 =5)   
         
plot(seqn1[,1:2],type="l",lty=1,xlab="Time",ylab="Prey");title("Stable focal point")
writelabel("A")
plot(seqn2[,1:2],type="l",xlab="Time",ylab="Prey")      ;title("Neutral stability")
writelabel("B")
plot(seqn3[,1:2],type="l",xlab="Time",ylab="species1")  ;title("Stable equilibrium")
writelabel("C")
plot(seqn4[,1:2],type="l",xlab="Time",ylab="species1")  ;title("Saddle point")
writelabel("D")

subtitle()


##############################
## equilibrium types        ##
##############################

# position of titles, number of figures, margins
par(mfrow=c(3,2),mar=c(3,4,3,1))
xlim          <- c(-1,1)            # axes limits
ylim          <- c(-1,1)

############################################
# Simple equation
############################################

equation <- function(ini,i1=5,a=-0.1,b=-0.3,cc=0,dd=0,verbose=FALSE,...)
{
eqn <- function (t,state,pars)
 {
  with (as.list(state),
  {
  dx<-a*x + cc*y
  dy<-b*y + dd*x
  list(c(dx,dy))
  })
 }
times <- seq(0,100,1)
state <- c(x=ini[1],y=ini[2])
out   <- as.data.frame(ode(state,times,eqn,parms=0))
lines(out$x,out$y)

Arrows(out$x[i1],out$y[i1],out$x[i1+1],out$y[i1+1],...)
#jacobian matrix: coefficients a,d,c,b
if (verbose) print(eigen(matrix(nrow=2,data=c(a,dd,cc,b)))$values)
 }

#---------------------------------------------#
# stable equilibrium
#---------------------------------------------#

plot(0,type="n",xlim=xlim,ylim=ylim,xlab="",ylab="",axes=FALSE,frame.plot=TRUE)
# Initial points outlined on a circle
x  <- seq(0,2*pi,by=pi/8)
xy <- cbind(1*cos(x), 1*sin(x))
for ( i in 1:nrow(xy)) equation(xy[i,],verbose=i==1)
points(0,0,pch=21,cex=3,bg="black",col="black")
title("Stable equilibrium")
writelabel("A")
#---------------------------------------------#
# unstable equilibrium
#---------------------------------------------#

plot(0,type="n",xlim=xlim,ylim=ylim,xlab="",ylab="",axes=FALSE,frame.plot=TRUE)
equation(c(0.01,-0.05),i1=20,a=0.2,b=0.2,cc=0.0,dd=0.2,verbose=TRUE)
equation(c(0.015,-0.05),i1=18,a=0.2,b=0.2,cc=0.0,dd=0.2)
equation(c(0.0,-0.05),i1=13,a=0.2,b=0.2,cc=0.0,dd=0.2)
equation(c(-0.01,0.05),i1=20,a=0.2,b=0.2,cc=0.0,dd=0.2)
equation(c(-0.015,0.05),i1=18,a=0.2,b=0.2,cc=0.0,dd=0.2)
equation(c(0.0,0.05),i1=13,a=0.2,b=0.2,cc=0.0,dd=0.2)

points(0,0,pch=21,cex=3,bg="white",col="black")
title("Unstable equilibrium")
writelabel("B")

#---------------------------------------------#
# saddle point
#---------------------------------------------#

plot(0,type="n",xlim=xlim,ylim=ylim,xlab="",ylab="",axes=FALSE,frame.plot=TRUE)
xy <- matrix(ncol=2,byrow=TRUE,
      data=c(-0.9,0.4,0.9,0.4,-0.9,-0.4,0.9,-0.4,
             -1,0,1,0,-0.8,0.2,0.8,0.2,-0.8,-0.2,0.8,-0.2))
for ( i in 1:nrow(xy)) equation(xy[i,],a=-0.1,b=0.1,verbose=i==1)
equation(c(0,0.05),a=-0.1,b=0.1,i1=25);equation(c(0,-0.05),a=-0.1,b=0.1,i1=25)
points(0,0,pch=21,cex=3,bg="darkgrey",col="black")
title("Saddle point")
writelabel("C")

#---------------------------------------------#
# limit cycles
#---------------------------------------------#

plot(0,type="n",xlim=xlim,ylim=ylim,xlab="",ylab="",axes=FALSE,frame.plot=TRUE)
xy <- matrix(ncol=2,byrow=TRUE,
      data=c(0,0.2,0,0.4,0,0.6,0,0.8,0,1))
for ( i in 1:nrow(xy)) equation(xy[i,],a=0,b=0,cc=-0.1,dd=0.1,verbose=i==1)
points(0,0,pch=21,cex=3,bg="lightgrey",col="lightgrey")
title("Neutral stability")
writelabel("D")

#---------------------------------------------#
# stable focal point
#---------------------------------------------#

plot(0,type="n",xlim=xlim,ylim=ylim,xlab="",ylab="",axes=FALSE,frame.plot=TRUE)
xy <- matrix(ncol=2,byrow=TRUE,
      data=c(0,0.2,0,0.4,0,0.6,0,0.8,0,1))
for ( i in 1:nrow(xy)) equation(xy[i,],a=0,b=-0.1,cc=-0.1,dd=0.1,verbose=i==1)
xy[,2]<- -xy[,2]
for ( i in 1:nrow(xy)) equation(xy[i,],a=0,b=-0.1,cc=-0.1,dd=0.1)
points(0,0,pch=21,cex=3,bg="black",col="black")

title("Stable focal point")
writelabel("E")

#---------------------------------------------#
# spirals
#---------------------------------------------#

plot(0,type="n",xlim=xlim,ylim=ylim,xlab="",ylab="",axes=FALSE,frame.plot=TRUE)
equation(c(0.05,-0.05),i1=40,a=0.,b=0.1,cc=0.1,dd=-0.1,verbose=TRUE)
equation(c(-0.05,0.05),i1=40,a=0.,b=0.1,cc=0.1,dd=-0.1)
equation(c(0.02,0.06),i1=40,a=0.,b=0.1,cc=0.1,dd=-0.1)
equation(c(-0.02,-0.06),i1=40,a=0.,b=0.1,cc=0.1,dd=-0.1)

points(0,0,pch=21,cex=3,bg="white",col="black")
title("Unstable focal point")
writelabel("F")
subtitle()

#############################################
## Fig. 7.10. types of limit cycles
#############################################

par(mfrow=c(2,2))
eqn <- function (t,state,pars)
 {
  with (as.list(c(state,pars)),
  {
  dx<-  a*y   +e*x*(x^2+y^2-1)
  dy<-  b*x   +f*y*(x^2+y^2-1)
  list(c(dx,dy))
  })
 }

equation2 <- function(ini,i1=5,a=-1,b=1,e=-1,f=-1,verbose=FALSE,endt=100,dt=0.1,...)
{
times <- seq(0,endt,dt)
state <- c(x=ini[1],y=ini[2])
out   <- as.data.frame(vode(state,times,eqn,parms=c(a=a,b=b,e=e,f=f)))
lines(out$x,out$y)
Arrows(out$x[i1],out$y[i1],out$x[i1+1],out$y[i1+1],...)
if(verbose)
{
Jacob    <- matrix(nrow=2,ncol=2,0)
equi <- c(x=0,y=0)
small    <- 1e-8
# reference model solution , at equilibirum values
ref  <- unlist(eqn(0,equi           ,pars=c(a=a,b=b,e=e,f=f)))

# increase sp1 and sp2 by very small amount
  pert1<- unlist(eqn(0,equi+c(small,0),pars=c(a=a,b=b,e=e,f=f)))
  pert2<- unlist(eqn(0,equi+c(0,small),pars=c(a=a,b=b,e=e,f=f)))

# finite differences
  Jacob[,1]<- (pert1-ref)/small
  Jacob[,2]<- (pert2-ref)/small

 # eigenvalues
  ei   <- eigen(Jacob)$values
print(ei)
}
return(out[nrow(out),2:3])
}

xlim<- c(-1.5,1.5)
ylim<- xlim

plot(0,type="n",xlim=xlim,ylim=ylim,xlab="",ylab="",axes=FALSE,frame.plot=TRUE)
equation2(c(-0.01,-0.01),i1=45,verbose=TRUE)
equation2(c(0.01,0.01),i1=45)
equation2(c(1.1,1.1),i1=5)
equation2(c(-1.1,1.1),i1=5)
equation2(c(-1.1,-1.1),i1=5)
equation2(c(1.1,-1.1),i1=5)
points(0,0,pch=21,cex=3,bg="lightgrey",col="lightgrey")
title("Stable limit cycle")


writelabel("A")



xlim<- c(-1.5,1.5)
ylim<- xlim

plot(0,type="n",xlim=xlim,ylim=ylim,xlab="",ylab="",axes=FALSE,frame.plot=TRUE)
equation2(c(-0.65,-0.65),i1=45,verbose=TRUE,endt=10,dt=0.01,e=1,f=1)
equation2(c(0.65,0.65),i1=45,endt=10,dt=0.01,e=1,f=1)
equation2(c(0.75,0.75),i1=45,endt=10,dt=0.01,e=1,f=1)
equation2(c(-0.75,-0.75),i1=45,endt=10,dt=0.01,e=1,f=1)
equation2(c(0.65,-0.65),i1=45,endt=10,dt=0.01,e=1,f=1)
equation2(c(0.75,-0.75),i1=45,endt=10,dt=0.01,e=1,f=1)
equation2(c(-0.65,0.65),i1=45,endt=10,dt=0.01,e=1,f=1)
equation2(c(-.75,0.75),i1=45,endt=10,dt=0.01,e=1,f=1)

points(0,0,pch=21,cex=3,bg="black",col="black")
title("Unstable limit cycle")

equation2(c(-0.71,0.70420315),i1=45,endt=8,dt=0.01,a=-1,b=1,e=1,f=1)


writelabel("B")
subtitle()


###############################################################################
####======================================================================#####
####                            R case studies                            #####
####======================================================================#####
###############################################################################

## The budworm model        ##
## Figs. 7.11, 7.12         ##

# parameter values
par(mfrow=c(2,2))
r   <- 0.05
K   <- 10    
bet <- 0.1   
alf <- 1     

# the model : density-dependent growth and sigmoid-type mortality rate
rate <- function(B,r=0.05) r*B*(1-B/K)-bet*B^2/(B^2+alf^2)

# plot the function for different values of r
Bseq <- seq(0,10,length=500)
rseq <- seq(0.01,0.07,by=0.02)
mat  <- outer (X=Bseq,Y=rseq, function(X,Y) rate(B=X,r=Y))
matplot(Bseq,mat,xlab="B",ylab="dB/dt",type="l",lty=1:10,col=1) 
abline(h=0,lty=2)                # add 0-axis
legend("bottomleft",legend=rseq,title="r",col=1,lty=1:10)
writelabel("A")

# equilibrium points, estimated numerically
equilibrium <- function(from=0,to=10,by=0.01,r)
{
Bseq <- seq(from,to,by=by)
mod  <- rate(B=Bseq,r=r)
Equi <- Bseq[which(mod==0)]

len  <- length(mod)
ss   <- mod[1:(len-1)]*mod[2:len]  # interval where functionvalues change sign
ii   <- which(ss<0)

root <- function(x) rate(B=x,r=r)
for (i in ii) Equi <- c(Equi,uniroot(root,lower=Bseq[i],upper=Bseq[i+1])$root)

return(Equi)
}

# Stability of equilibrium points
stability <- function(equi,r)
{
tiny      <- 1e-8
ref       <- rate(equi,r=r)
eigenval  <-sign(((rate(B=equi+tiny,r=r)-ref)/tiny) ) # sign of eigenvalue
return(eigenval)
}

curve(rate(x,0.05),ylab="dB/dt",main="r=0.05",from=0,to=10)
abline(h=0)
Eq    <- equilibrium(r=0.05)
eig   <- stability(Eq,0.05)
points(x=Eq,y=rep(0,length(Eq)),pch=21,cex=2,
bg=c("grey","black","white")[eig+2] )
writelabel("B")


# bifurcation diagram - fig. 7.12
rseq <- seq(0.01,0.07,by=0.0001)

plot(0,xlim=range(rseq),ylim=c(0,10),type="n",
xlab="r",ylab="B*",main="spruce budworm model")

for (r in rseq) {
equi <- equilibrium(r=r)
eig  <- stability(equi,r) 

  points(rep(r,length(equi)),equi,pch=22,
   col=c("darkgrey","black","lightgrey")[eig+2],
   bg =c("darkgrey","black","lightgrey")[eig+2]) 
}

legend("topleft",pch=22,pt.cex=2,c("stable","unstable"),
col=c("black","lightgrey"),pt.bg=c("black","lightgrey"))

equi <-equilibrium(r=0.05) 
arrows(0.05,10         ,0.05,equi[4]+0.2,length=0.1 )
arrows(0.05,equi[3]+0.2,0.05,equi[4]-0.2,length=0.1 )
arrows(0.05,equi[3]-0.2,0.05,equi[2]+0.2,length=0.1 )
arrows(0.05,equi[1]+0.2,0.05,equi[2]-0.2,length=0.1 )  

equi <-equilibrium(r=0.038) 
arrows(0.038,10         ,0.038,equi[2]+0.1,length=0.1 )
arrows(0.038,equi[1]+0.1,0.038,equi[2]-0.1,length=0.1 )

equi <-equilibrium(r=0.07) 
arrows(0.07,10         ,0.07,equi[2]+0.2,length=0.1 )
arrows(0.07,equi[1]+0.2,0.07,equi[2]-0.2,length=0.1 )
subtitle()

par(mfrow=c(1,1))

##############################
## Case study: competition  ##
##############################

r1    <- 3              # parameters
r2    <- 2
K1    <- 1.5
K2    <- 2
alf12 <-1
alf21 <-2

Lotka<-function(t,N,pars)

{
 dN1 <- r1*N[1]*(1-(N[1]+alf12* N[2])/K1)
 dN2 <- r2*N[2]*(1-(N[2]+alf21* N[1])/K2)

 list(c(dN1  , dN2  ))     # the rate of change
}

Ax  <- c(0,K2/alf21)
Ay  <- K2 - alf21* Ax
By  <- c(0,K1/alf12)
Bx  <- K1 - alf12* By
xlim   <- range(c(Ax, Bx))
ylim   <- range(c(Ay, By))

plot  (x=Ax,y=Ay, type="l", lwd=3,   # 1st isocline
     main="Competition phase-plane",
       xlab="N1",ylab="N2",xlim=xlim,ylim=ylim)
lines (Bx,By,lwd=3,lty=2)            # 2nd isocline

tarrows <- function(out,ds,...)
{
 # select point at predefined distance from point (p1,p2)
 p   <- unlist(out[nrow(out),2:3])
 dd  <- (out[,2]-p[1])^2+ (out[,3]-p[2])^2
 dd2 <-  c(dd[-1],dd[1])
 i1<-which(dd<ds&dd2>ds | dd>ds&dd2<ds)

 p   <- unlist(out[1,2:3])
 dd  <- (out[,2]-p[1])^2+ (out[,3]-p[2])^2
 dd2 <-  c(dd[-1],dd[1])
 i2<-which(dd<ds&dd2>ds | dd>ds&dd2<ds)[1]
 ii <- c(i1,i2)
 for (i in ii ) arrows(out[i,2],out[i,3],out[i+1,2],out[i+1,3],length=0.1,lwd=1,...)

}

trajectory <- function(N1,N2)
{
times  <-seq(0,30,0.1)
state  <-c(N1 = N1, N2 = N2)
out    <-as.data.frame(ode(state,times,Lotka,0))

lines (out$N1,out$N2,type="l")
tarrows(out,ds)
}

ds<- 0.1

trajectory (0.05,0.3)
trajectory (0.11,0.3)
trajectory (1.5,1.8)
trajectory (1.0,2.0)

# 4 equilibrium points
X  <- c(0,0 ,K1,(K1-alf12*K2)/(1-alf12*alf21))
Y  <- c(0,K2,0 ,(K2-alf21*K1)/(1-alf12*alf21))

# Jacobian matrix, and eigenvalues
small <- 1e-8
Jacob <- matrix(nrow=2,ncol=2)
ei    <- matrix(nrow=4,ncol=2)

for ( i in 1:4)
{
 N1 <- X[i]
 N2 <- Y[i]
 # reference model solution
 ref  <- unlist(Lotka(0,c(N1=N1,N2=N2)      ,0))

 # increase sp1 and sp2 by very small amount
 pert1<- unlist(Lotka(0,c(N1=N1+small,N2=N2),0))
 pert2<- unlist(Lotka(0,c(N1=N1,N2=N2+small),0))

 # finite differences
 Jacob[,1]<- (pert1-ref)/small
 Jacob[,2]<- (pert2-ref)/small

 # eigenvalues
 ei[i,]   <- eigen(Jacob)$values

 # white:unstable node, black:stable node, grey:saddle

 if (sign(ei[i,1])>0 & sign(ei[i,2])>=0) col <- "white"
 if (sign(ei[i,1])<0 & sign(ei[i,2])<=0) col <- "black"
 if (sign(ei[i,1])* sign(ei[i,2])   <0 ) col <- "grey"

# equilibrium point plotting
 points(N1,N2,pch=22,cex=2.0,bg=col,col="black")

}
cbind(N1=X,N2=Y,ei)

eig      <- eigen(Jacob)
vv       <- eig$vector[eig$values<0]



# the reverse model, output -rate of change
revmod <- function(t,N,p) list(-1*unlist(Lotka(t,N,p)))
times  <-seq(0,1.9,0.05)

# first direction
state <- c(N1,N2) + 0.01*vv

out   <-as.data.frame(ode(state,times,revmod,0))
lines(out[,2],out[,3],lty=2)
tarrows(out,ds,code=1)

# second direction
state <- c(N1,N2) - 0.01*vv
times  <-seq(0,10,0.05)
out   <-as.data.frame(ode(state,times,revmod,0))
lines(out[,2],out[,3],lty=2)
tarrows(out,ds,code=1)

# trajectories out of the equilibrium point
ww       <- eig$vector[eig$values>0]
trajectory (N1+0.05*ww[1],N2+0.05*ww[2])
trajectory (N1-0.05*ww[1],N2-0.05*ww[2])

legend("right",legend=c("isocline N1","isocline N2","trajectory","separatrice",
"saddle point","stable equilibrium","unstable equilibrium"),lty=c(2,1,1,2,NA,NA,NA),
lwd=c(2,2,1,1,NA,NA,NA),pch= c(NA,NA,NA,NA,22,22,22),pt.bg=c(NA,NA,NA,NA,"grey","black","white"))

subtitle()


##############################
## The Lorenz equations     ##
## chaos                    ##
##############################

#----------------------#
# the model equations: #
#----------------------#

Lorenz<-function(t,state,parameters)
  {
  with(as.list(c(state)),{ 

    ## the rates of change of state variables
    dx     <- -8/3*x+y*z
    dy     <- -10*(y-z)
    dz     <- -x*y+28*y-z

    ## the output
    list(c(dx,dy,dz))            })
 }  # end of model

#-------------------------#
# the initial conditions: #
#-------------------------#
 
state     <-c(x=1,
              y=1,
              z=1)

#----------------------#
# RUNNING the model:   #
#----------------------#

times <-seq(0,100,0.001)
out   <-as.data.frame(ode(state,times,Lorenz,0))

#------------------------#
# PLOTTING model output: #
#------------------------#
par(mfrow=c(1,1))
require(scatterplot3d)
scatterplot3d(out$x,out$y,out$z,type="l",main="Lorenz butterfly",ylab="",
              grid=FALSE,box=FALSE)

subtitle()



###############################################################################
# Silicate diagenesis
###############################################################################
require(rootSolve)
SiDiamodel <- function (time=0,Conc,pars=NULL)
{
 BSi<- Conc[1:N]
 DSi<- Conc[(N+1):(2*N)]

######################
# transport          #
######################
# diffusive fluxes at upper interface of each layer

# upper concentration imposed (bwDSi), lower: zero gradient
 DSiFlux <- -SedDisp *   IntPor *diff(c(bwDSi ,DSi,DSi[N]))/thick    
 BSiFlux <- -Db      *(1-IntPor)*diff(c(BSi[1],BSi,BSi[N]))/thick 

 BSiFlux[1] <- BSidepo                         # upper boundary flux is imposed

######################
# BSi dissolution    #
######################
 Dissolution <- rDissSi * BSi*(1.- DSi/EquilSi )^pow 
 Dissolution <- pmax(0,Dissolution)

# Rate of change= Flux gradient, corrected for porosity and dissolution
 dDSi        <- -diff(DSiFlux)/thick/Porosity      +           # transport
                 Dissolution * (1-Porosity)/Porosity           # biogeochemistry

 dBSi        <- -diff(BSiFlux)/thick/(1-Porosity)  - Dissolution				

 return(list(c(dBSi=dBSi,dDSi=dDSi),
             Dissolution=Dissolution,
             DSiSurfFlux =DSiFlux[1],DSIDeepFlux =DSiFlux[N+1],
             BSiDeepFlux =BSiFlux[N+1]))
}

# sediment parameters
thick    <- 0.05                       # thickness of sediment layers (cm)
Intdepth <- seq(0,10,by=thick)         # depth at upper interface of each layer 
Nint     <- length(Intdepth)           # number of interfaces
Depth    <- 0.5*(Intdepth[-Nint] +Intdepth[-1]) # depth at middle of each layer
N        <- length(Depth)                       # number of layers

por0    <- 0.9                         # surface porosity (-)
pordeep <- 0.7                         # deep porosity    (-)
porcoef <- 2                           # porosity decay coefficient  (/cm)
Porosity <- pordeep + (por0-pordeep)*exp(-Depth*porcoef)     # porosity profile, middle of layers
IntPor   <- pordeep + (por0-pordeep)*exp(-Intdepth*porcoef)  # porosity profile, upper interface

dB0      <- 1/365           # cm2/day       - bioturbation coefficient
dBcoeff  <- 2               
mixdepth <- 5                # cm
Db       <- pmin(dB0,dB0*exp(-(Intdepth-mixdepth)*dBcoeff))

# biogeochemical parameters
SedDisp  <- 0.4              # diffusion coefficient, cm2/d  
rDissSi  <- 0.005            # dissolution rate, /day
EquilSi  <- 800             # equilibrium concentration
pow      <- 1
BSidepo  <- 0.2*100          # nmol/cm2/day
bwDSi    <- 150              # mumol/l  

# initial guess of state variables-just random numbers between 0,1
Conc     <- runif(2*N) 

# three runs with different deposition rates
BSidepo  <- 0.2*100          # nmol/cm2/day
sol  <- steady.1D (y=Conc, func=SiDiamodel, parms=NULL, nspec=2) 
CONC <- sol$y

BSidepo  <- 2*100          # nmol/cm2/day
sol2 <- steady.1D (y=Conc, func=SiDiamodel,parms=NULL, nspec=2) 
CONC2 <- sol2$y

BSidepo  <- 3*100          # nmol/cm2/day
sol3 <- steady.1D (y=Conc, func=SiDiamodel,parms=NULL, nspec=2) 
CONC3 <- sol3$y

DSi  <- cbind(CONC[(N+1):(2*N)],CONC2[(N+1):(2*N)],CONC3[(N+1):(2*N)])
BSi  <- cbind(CONC[1:N],CONC2[1:N],CONC3[1:N])

par(mfrow=c(2,2))

matplot(DSi,Depth,ylim=c(10,0),xlab="mmolSi/m3 Liquid",main="DSi",type="l",lwd=c(1,2,1),col="black")
matplot(BSi,Depth,ylim=c(10,0),xlab="mmolSi/m3 Solid" ,main="BSi",type="l",lwd=c(1,2,1),col="black")
legend("right",c("0.2","2","3"),title="mmol/m2/d",lwd=c(1,2,1),lty=1:3)
plot(Porosity,Depth,ylim=c(10,0),xlab="-" ,main="Porosity",type="l",lwd=2)
plot(Db,Intdepth,ylim=c(10,0),xlab="cm2/d" ,main="Bioturbation",type="l",lwd=2)



##############################
## Marine zooplankton in    ##
## estuary, 1st order decay ##
##############################


Estuary <-function(t,C,pars)
{

  with (as.list(pars),{

    Input <-  Q     * c(Criver, C) +
             -Estar * diff(c(Criver, C, Csea))
    dC    <- -diff(Input)/Volume + rate*C
    list(dC)
  })
}


nbox    <- 100
Length  <- 100000                       # m total length of estuary
dx      <- Length/nbox

# distance from river to interfaces and to mid of boxes  (m)
IntDist   <- seq(0,by=dx,length.out=nbox+1)
Dist      <- seq(dx/2,by=dx,length.out=nbox)

# Sigmoidal increase in cross sectional area
# Cross section at interfaces and midle of boxes         (m2)
Area_int <- 4000 + 76000 * IntDist^5 /(IntDist^5+50000^5)
Area     <- 4000 + 76000 * Dist^5    /(Dist^5+50000^5)

# Volume of boxes                                        (m3)
Volume   <- Area*dx

Eriver   <- 0                      # m2/d dispersion coefficient at upstream boundary
Esea     <- 350*3600*24            # m2/d dispersion coefficient at seaward boundary
E        <- Eriver + IntDist/Length * Esea

# BulkDisp is used in the calculations; m3/d
Estar  <- E * Area_int/dx

# the model parameters:
pars   <-   c(Criver  = 0.0,          # riverine conc
              Csea    = 100.0,        # seaward conc
              rate    = -0.05,        # /day  growth rate
              Q       = 100*3600*24)  # m3/d, mean river flow

pars["rate"] <- 0
pars["Csea"] <- 35

sal <- steady.band(y=rep(35,times=nbox),func=Estuary,
                   parms=pars,nspec=1)$y

Zoo  <- NULL

gSeq <-seq (-0.05,0.01, by=0.01)
for (g in gSeq)
{
pars["rate"] <- g
pars["Csea"] <- 100
st <- steady.band(y=rep(100,nbox),func=Estuary,parms=pars,nspec=1 ,pos=TRUE  ) #
Zoo <-  cbind(Zoo,st$y)
}
par(mfrow=c(1,1))
matplot(sal,Zoo,type="l",lwd=1,col="black",xlab="Salinity",
        ylab="gDWT/m3",main="Zooplankton",lty=1:20)
legend("topleft",legend= gSeq,title="g, /day",lty=1:20)
subtitle()

###############################################################################
####======================================================================#####
####                               Projects                               #####
####======================================================================#####
###############################################################################
## Schaefer fisheries model ##
##############################

par(mfrow=c(2,2))
# parameter values
r   <- 0.05
K   <- 10    
mrt <- 0.02   

# the model : density-dependent growth and sigmoid-type mortality rate
rate <- function(x,mrt=0.02) r*x*(1-x/K)-mrt*x

# plot the function for different values of mortality
xseq   <- seq(0,10,length=500)
mrtseq <- seq(0.0,0.1,by=0.02)
mat  <- outer (X=xseq,Y=mrtseq, function(X,Y) rate(x=X,mrt=Y))
matplot(xseq,mat,xlab="x",ylab="dx/dt",type="l",lty=1:10,col=1,main="Schaefer model") 
abline(h=0,lty=2)                # add 0-axis
legend("bottomleft",legend=mrtseq,title="mrt",col=1,lty=1:10)
writelabel("A")


par(mar=c(5.1,4.1,4.1,4.1))
mrtseq <- seq(0.0,0.05,by=0.0005)
equi  <- K*(1-mrtseq/r)
Yield <- equi*mrtseq
plot(mrtseq,equi,type="l",lwd=2,xlab="fishing mortality",ylab="Equilibrium biomass", 
     main="Schaefer model")
par(new=TRUE)
plot(mrtseq,Yield,axes=FALSE,xlab="",ylab="",type="l",lwd=1)
axis(4)
mtext(side=4,line=2.5,"Fisheries yield",cex=0.8)
mmax <-mrtseq[which.max(Yield)]
xmax <- equi  [which.max(Yield)]

c(mmax,xmax)
 
abline(v=mmax,lty=3)
legend("topright",lwd=c(1,2),legend=c("Yield","Biomass"))


r<-2.61
K <- 1.34e8
Q<-3.8e-5
Estar <- 0.5*r/Q
Estar

writelabel("B")

par(mar=c(5.1,4.1,4.1,2.1))

## Fisheries model with     ##
## Allee effect             ##
##############################

# parameter values
r   <- 0.05
K   <- 10    
K0  <- 1.

# the model : density-dependent growth, allee effect and linear mortality
rate <- function(x,mrt=0.02) r*x*(1-x/K)*(x/K0-1)-mrt*x

# plot the function for different values of mortality
xseq   <- seq(0,10,length=500)
mrtseq <- seq(0.0,0.1,by=0.02)
mat  <- outer (X=xseq,Y=mrtseq, function(X,Y) rate(x=X,mrt=Y))
matplot(xseq,mat,xlab="x",ylab="dx/dt",type="l",lty=1:10,col=1,main="fisheries model with Allee effect") 
abline(h=0,lty=2)                # add 0-axis
legend("bottomleft",legend=mrtseq,title="mrt",col=1,lty=1:10)
writelabel("C")

# equilibrium points by numerical approximation
equilibrium <- function(from=0,to=10,by=0.01,mrt)
{
xseq <- seq(from,to,by=by)
mod  <- rate(x=xseq,mrt=mrt)
Equi <- xseq[which(mod==0)]

len  <- length(mod)
ss   <- mod[1:(len-1)]*mod[2:len]  # interval where functionvalues change sign
ii   <- which(ss<0)

root <- function(x) rate(x=x,mrt=mrt)
for (i in ii) Equi <- c(Equi,uniroot(root,lower=xseq[i],upper=xseq[i+1])$root)

return(Equi)
}

# Stability of equilibrium points
stability <- function(equi,mrt)
{
tiny      <- 1e-8
ref       <- rate(x=equi,mrt=mrt)
eigenval  <-sign(((rate(x=equi+tiny,mrt=mrt)-ref)/tiny) ) # sign of eigenvalue
return(eigenval)
}

# bifurcation diagram
mrtseq <- seq(0.0,0.15,by=0.0005)
plot(0,xlim=range(mrtseq),ylim=c(0,10),type="n",
xlab="fishing mortality",ylab="Equilibrium biomass", 
     main="fisheries model with Allee effect")

for (mrt in mrtseq) {
equi <- equilibrium(mrt=mrt)
eig  <- stability(equi,mrt) 

  points(rep(mrt,length(equi)),equi,pch=22,
   col=c("darkgrey","black","lightgrey")[eig+2],
   bg =c("darkgrey","black","lightgrey")[eig+2]) 
}

legend("topright",pch=22,pt.cex=2,c("stable","unstable"),
col=c("darkgrey","lightgrey"),pt.bg=c("darkgrey","lightgrey"))
writelabel("D")
subtitle()

##############################
## predator-prey with       ##
## type-II response         ##
##############################

########################################################################
# Stability properties of model with density dependent growth of prey
# and holling type II functional response of grazing
########################################################################

# parameters
r     <- 3      # prey rate of increase        
K     <- 10     # prey carrying capacity
ks    <- 2      # half-saturation constant
e     <- 0.5    # growth efficiency of predator
m     <- 0.2    # predator mortality rate
In    <- 0.6    # ingestion rate predator

#========================================
# rate of change 
#========================================
model<-function(t,State,pars)
{
 with (as.list(State),{
 dN <- r*Prey*(1-Prey/K)-In*Prey/(Prey+ks)*Predator
 dP <- e*In*Prey/(Prey+ks)*Predator - m*Predator

 list(c(dN  , dP  ))     # the rate of change
 })
}

#========================================
# isoclines and roots
#========================================
isocline <- function()
{

# prey isocline is a curve
 curve  (r/In*(1-x/K)*(x+ks), from=0, to= K,lwd=3,   # 1st isocline
        xlab="Prey",ylab="Predator",xlim=xlim,ylim=ylim)

# predator isocline is a vertical line
 Neq    <- m*ks/(In*e-m)        
 if (Neq < 0 | is.infinite(Neq)) {warnings("Model has no solution");return}
 abline (v=Neq,lwd=3,lty=2)                # 2nd isocline
# legend("topright",legend=paste("r=",In))

# first equilibrium point : (0,0)
# then intersection between predator isocline and X-axis, 
# and intersection between predator and prey isocline
 roots <- matrix(ncol=2,byrow=TRUE,data=
          c(0  ,0,
            K  ,0,
            Neq, r/In*(1-Neq/K)*(Neq+ks)))
 colnames(roots) <- c("Prey","Predator")

 root<-stability(roots)
 return(root)
 
}
# end isocline

#========================================
# the stability properties of the equilibrium points
#========================================

stability <- function (roots)
{
 small    <- 1e-8
 # Jacobian matrix, and eigenvalues
 Jacob  <- matrix(nrow=2,ncol=2,0)
 eig    <- NULL
 for (i in 1:nrow(roots))
 {

  equi <- roots[i,]

# reference model solution , at equilibirum values
  ref  <- unlist(model(0,equi           ,0))

# increase sp1 and sp2 by very small amount
  pert1<- unlist(model(0,equi+c(small,0),0))
  pert2<- unlist(model(0,equi+c(0,small),0))

# finite differences
  Jacob[,1]<- (pert1-ref)/small
  Jacob[,2]<- (pert2-ref)/small

 # eigenvalues
  ei   <- eigen(Jacob)$values
  eig <- rbind(eig,ei)

 # white:unstable node, black:stable node, grey:saddle
  ifelse (is.complex(ei), pch<-21, pch<-22)
  ei <- Re(ei)
  if (sign(ei[1])>0 & sign(ei[2])>=0) pcol <- "white"
  if (sign(ei[1])<0 & sign(ei[2])<=0) pcol <- "black"
  if (sign(ei[1])   * sign(ei[2])<0 ) pcol <- "grey"
 # equilibrium point plotting

  points (equi[1],equi[2],cex=2,pch=pch,bg=pcol,col="black")

 }
 return(list(equilibrium=roots,eigenvalues=eig) )
}   # end stability


#========================================
# model trajectories
#========================================

trajectory <- function(N,P)
{
times  <-seq(0,100,0.1)
state  <-c(Prey = N, Predator = P)
out    <-as.data.frame(ode(state,times,model,0))

lines (out$Prey,out$Predator,type="l")
arrows(out[10,2],out[10,3],out[11,2],out[11,3],length=0.1,lwd=1)

}

#========================================
# model applications
#========================================

par(oma=c(0,0,1,0),mfrow=c(2,2))

In     <- 0.55    # ingestion rate predator
xlim   <- c(0,K)
ylim   <- c(0,25)
isocline( ) 
trajectory (6,10)
trajectory (10,15)
trajectory (4,5)
trajectory (2,20)
title("stable equilibrium")

In     <- 0.45    # ingestion rate predator
ylim   <- c(0,30)
xlim   <- c(0,20)
isocline( ) 
trajectory (6,5)
trajectory (10,25)
trajectory (20,20)

title("stable equilibrium")

In     <- 0.8    # ingestion rate predator
ylim   <- c(0,20)
xlim   <- c(0,K)
isocline( ) 
trajectory (6,8)
trajectory (10,5)
trajectory (1,1)
trajectory (2.1,12.1)
title("stable limit cycle")

mtext(3,text="Predator-Prey, Monod grazing",line=-1,outer=TRUE,cex=1.2)


plot(0,type="n",axes=FALSE,xlab="",ylab="")
legend("center",pch=c(22,22,21,NA,NA),pt.bg=c("black","grey","white",NA,NA),
      lty=c(NA,NA,NA,1,2),lwd=c(NA,NA,NA,2,2),pt.cex=2,title="equilibrium",
      legend=c("stable","unstable","neutral","constant prey","constant predator"))
subtitle()

###############################################################################
# NPZmodel equations
###############################################################################


NPZmodel <- function (time=0,state,pars=NULL)
{
 N<- state[1       :Nb    ]
 P<- state[(Nb+1)  :(2*Nb)]
 Z<- state[(2*Nb+1):(3*Nb)]

# transport          #
# advective fluxes at upper interface of each layer
# freshwater concentration imposed 

 NFlux <- flow * c(RiverN ,N) 
 PFlux <- flow * c(RiverP ,P) 
 ZFlux <- flow * c(RiverZ ,Z) 
 
# Biology            #

 Pprod <-  mumax * N/(N+kN)*P      # primary production
 Graz  <-  gmax * P/(P+kP)*Z       # zooplankton grazing
 Zmort <-  mrt*Z                   # zooplankton mortality

# Rate of change= Flux gradient and biology
 dN        <- -diff(NFlux)/delx - Pprod + Graz*(1-eff) + Zmort
 dP        <- -diff(PFlux)/delx + Pprod - Graz 
 dZ        <- -diff(ZFlux)/delx         + Graz*eff     - Zmort 

 return(list(c(dN=dN,dP=dP,dZ=dZ),
             c(Pprod=Pprod,Graz=Graz,Zmort=Zmort,
            Nefflux=NFlux[Nb],Pefflux=PFlux[Nb],Zefflux=ZFlux[Nb])))
}


###############################################################################
# model application
###############################################################################
# model parameters
riverlen <- 100               # total length river,         km
Nb       <- 100               # number boxes
delx     <- riverlen / Nb     # box length

RiverN   <- 100               # N concentration at river,   mmolN/m3
RiverP   <- 10                # P concentration at river,   mmolN/m3
RiverZ   <- 1                 # Z concentration at river,   mmolN/m3

flow     <- 1                  # river flow,                 km/day

mumax    <- 0.5               # maximal light-limited primary production, /day
kN       <- 1                 # half-saturated N for pprod, mmolN/m3
gmax     <- 0.5               # max grazing rate,           /day
kP       <- 1                 # half-saturated P for graz , mmolN/m3
mrt      <- 0.05               # mortality rate,             /day
eff      <- 0.7               # growth efficiency,          -

# initial guess of NPZ; just random numbers
Conc     <- runif(3*Nb) #c(N,P,Z)

# steady-state solution, v= 1 km/day
sol  <- steady.1D (y=Conc, func=NPZmodel,nspec=3,maxiter=100,atol=1e-10,positive=TRUE) 
Conc <- sol$y
N    <- Conc[1       :Nb    ]
P    <- Conc[(Nb+1)  :(2*Nb)]
Z    <- Conc[(2*Nb+1):(3*Nb)]

Res  <- NPZmodel(state=Conc)

# run with v=5 km d-1
flow     <- 5                  # river flow,                 km/day
Conc     <- runif(3*Nb) #c(N,P,Z)
sol2  <- steady.1D (y=Conc, func=NPZmodel,nspec=3,maxiter=100,atol=1e-10,positive=TRUE) 
Conc2 <- sol2$y
N2    <- Conc2[1       :Nb    ]
P2    <- Conc2[(Nb+1)  :(2*Nb)]
Z2    <- Conc2[(2*Nb+1):(3*Nb)]

# run with v=10 km d-1
flow     <- 10                  # river flow,                 km/day
Conc     <- runif(3*Nb) #c(N,P,Z)
sol3  <- steady.1D (y=Conc, func=NPZmodel,nspec=3,maxiter=100,atol=1e-10,positive=TRUE) 
Conc3 <- sol3$y
N3    <- Conc3[1       :Nb    ]
P3    <- Conc3[(Nb+1)  :(2*Nb)]
Z3    <- Conc3[(2*Nb+1):(3*Nb)]


par(mfrow=c(2,2))
Dist <- seq(delx/2,riverlen,by=delx)
plot(Dist,N,ylab="mmolN/m3",main="DIN",type="l",lwd=2,ylim=c(0,110))
lines(Dist,N2)
lines(Dist,N3,lty=2)
plot(Dist,P,ylab="mmolN/m3" ,main="Phytoplankton",type="l",lwd=2,ylim=c(0,110))
lines(Dist,P2)
lines(Dist,P3,lty=2)
plot(Dist,Z,ylab="mmolN/m3" ,main="Zooplankton",type="l",lwd=2,ylim=c(0,110))
lines(Dist,Z2)
lines(Dist,Z3,lty=2)

plot(0,type="n",xlab="",ylab="",axes=FALSE)
legend("center",lty=c(1,1,2),lwd=c(2,1,1),
c(expression(v==1~km~d^{-1}),expression(v==5~km~d^{-1}),expression(v==10~km~d^{-1}))) 

mtext(side=3,outer=TRUE,"NPZ model",line=-1.5,cex=1.5)
subtitle()

par(ask=opar$ask)
par(mar=opar$mar)
par(oma=opar$oma)
par(mfrow=opar$mfrow)

