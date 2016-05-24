##############################
## Soetaert and Herman      ##
## ecological modelling     ##
## Figures from chapter 10  ##
## Dynamic programming      ##
##############################

opar <- par()
par(ask=TRUE)
par(mfrow=c(1,1))
figNr <- 1
subtitle <- function()
{
 mtext(side=1,outer=TRUE,"Soetaert and Herman - chapter 10  ",cex=0.7,adj=1,line=-1.5)
 mtext(side=1,outer=TRUE,paste("  Fig. 10.",figNr,sep=""),cex=0.7,adj=0,line=-1.5)
 figNr <<- figNr +1
}

###############################################################################
####======================================================================#####
####                                Theory                                #####
####======================================================================#####
###############################################################################

########################################
# Figure 10.1. General dynamic programming scheme
########################################

par(mar=c(2,2,2,2))
openplotmat()
elpos <-coordinates (c(1,1,1,1,1,1))
segmentarrow(elpos[5,],elpos[2,],arr.side=3,arr.pos=0.15,dd=0.5)
dx  <- 0.35
dy  <- 0.02
for (i in c(1:5)) straightarrow (elpos[i,],elpos[i+1,],lwd=2,arr.pos=0.5)
textround(elpos[1,],dx,2*dy,lab=c("Final fitness for all states X","f(X,Tend)"))
textround(elpos[2,],dx,2*dy,lab=c("Fitness for all states X and choices i",
         "Vi(X,T) = direct fitness gain + expected gain of next state"))
textround(elpos[3,],dx,2*dy,lab=c("Optimal fitness for all states",
                               "f(X,T) = max Vi(X,T)") )
textround(elpos[4,],dx,dy,lab="T = T-1")
textdiamond(elpos[5,],0.07,lab="T=0?")
textrect(elpos[6,],dx,dy ,lab="Optimal policy")

text(0.2,elpos[5,2]+0.015,"No",adj=0,font=3)
text(0.52,elpos[5,2]-0.08,"yes",adj=0,font=3)
subtitle()



#################################################################
## Fig. 10.2 optimal selection of FEEDING OR REPRODUCING patch ##
#################################################################

################################
# parameters settings          #
################################

x_crit  <- 5                 # critical mass to survive
x_max   <- 15               # maximal mass
x_rep   <- 6                 # critical mass for reproduction
x_class <- x_crit:x_max      # biomass classes
nmass   <- length(x_class)   # number of mass classes

t_max   <- 20               # number of time steps
times   <- 1:(t_max-1)
npatch  <- 2                 # number of patches

cost     <- c(1, 1)          # cost of a patch
feedgain <- c(2 ,0)          # gain of feeding
repr     <- c(0 ,1)          # max reproduction
psurv    <- c(0.9,0.95)      # survival probability

################################
# result matrices              #
################################

f         <- matrix(nrow=t_max,ncol=nmass ,0)    # optimal fitness values
bestpatch <- matrix(nrow=t_max-1,ncol=nmass-1,0) # best patch choice
V         <- vector(length=npatch)               # current fitness for a patch

#####################################
# find fitness for mass x at time t #
#####################################

fitness <- function(x,t)

{
 xx <- pmin(x ,x_max)
 xx <- pmax(xx,x_crit)
 fitness <- f[t,xx-x_crit+1]
}


################################
# main optimisation loop       #
################################

optimizePatch <- function ()
{
for (t in rev(times))                         # backward in time
{
 for (x in x_class[-1])                       # for each biomass class, except x-crit
 {
  dfit <- pmax(0,pmin(x-x_rep,repr))          # reproduction
  expectgain <- psurv * fitness(x-cost+feedgain-dfit,t+1)
  V          <- dfit + expectgain

  f[t,x-x_crit+1]      <<- max(V)             # optimal fitness

  bestpatch[t,x-x_crit] <- which.max(V)[1]    # best patch

}                                             # next biomass class x
}                                             # next time t
 return(list(bestpatch=bestpatch, fitness=f))
} # end function optimizePatch

###############################################
# RUN 1: animal dies at end
###############################################

# final fitness, at t_max
fend      <- 0
f[t_max,] <- fend

Opt <- optimizePatch()

par(mar=c(5.1,4.1,4.1,2.1),mfrow=c(2,2))
image(x=times,y=x_class[-1],z=Opt$bestpatch,ylab="Biomass",xlab="time",
      main="optimal strategy - semelparous",col=c("darkgrey","white"),zlim=c(1,2))
box()
contour(x=1:t_max,y=x_class,z=Opt$fitness,add=TRUE,labcex=1.1)
legend("topright",fill=c("darkgrey","white"),legend=c("feeding","reproducing"),bg="white")

###############################################
# RUN 2: animal fitness at end ~ animal biomass
###############################################

# final fitness, at t_max
fend      <- 0:(nmass-1)
f[,]      <- 0
f[t_max,] <- fend

Opt2<- optimizePatch()

image(x=times,y=x_class[-1],z=Opt2$bestpatch,ylab="Biomass",xlab="time",
      main="optimal strategy - iteroparous",col=c("darkgrey","white"),zlim=c(1,2))
box()
contour(x=1:t_max,y=x_class,z=Opt2$fitness,add=TRUE,labcex=1.1)
legend("topright",fill=c("darkgrey","white"),legend=c("feeding","reproducing"),bg="white")

subtitle()


###############################################################################
####======================================================================#####
####                             R case study                             #####
####======================================================================#####
###############################################################################


#############################################################
## THE PATCH SELECTION MODEL                               ##
## Clark CW and M Mangel 1999.                             ##
## Dynamic State Variable Models in Ecology: methods and   ##
## applications. Oxford University Press.                  ##
## Chapter 1.                                              ##  
#############################################################


################################
# parameters settings          #
################################

x_crit  <- 0                 # critical mass to survive
x_max   <- 30                # maximal mass
x_rep   <- 4                 # critical mass for reproduction
x_class <- x_crit:x_max      # biomass classes
nmass   <- length(x_class)   # number of mass classes   

t_max   <- 20                # number of time steps 
times   <- 1:(t_max-1)
npatch  <- 3                 # number of patches

psurvive <- c(0.99,0.95,0.98)   # probability of surviving
pfood    <- c(0.2 ,0.5 ,0 )     # probability of feeding
cost     <- c(1 ,1 ,1 )         # cost of a patch
feedgain <- c(2 ,4 ,0 )         # gain of feeding
repr     <- c(0 ,0 ,4 )         # max reproduction

################################
# result matrices              #
################################

f         <- matrix(nrow=t_max,ncol=nmass ,0)    # optimal fitness values
bestpatch <- matrix(nrow=t_max-1,ncol=nmass-1,0) # best patch choice
V         <- vector(length=npatch)               # current fitness for a patch

# final fitness, at t_max
fend      <- 60
kx        <- 0.25*x_max

f[t_max,]  <- fend*(x_class-x_crit)/(x_class-x_crit+kx)

#####################################
# find fitness for mass x at time t #
#####################################

fitness <- function(x,t)

{
xx <- pmin(x ,x_max) 
xx <- pmax(xx,x_crit) 
fitness <- f[t,xx-x_crit+1]
}


################################
# main optimisation loop       #
################################


for (t in rev(times))               # backward in time
{
 for (x in x_class[-1])             # for each biomass class, except x-crit
 {                                  
  dfit <- pmax(0,pmin(x-x_rep,repr)) # reproduction
  expectgain <- psurvive*( pfood   *fitness(x-cost+feedgain-dfit,t+1) + 
                          (1-pfood)*fitness(x-cost-dfit ,t+1) )
  V    <- dfit + expectgain  
  V[expectgain == 0] <- 0           # dead
  f[t,x-x_crit+1]       <- max(V)          # optimal fitness
  bestpatch[t,x-x_crit] <- which.max(V)    # best patch 
       
 }                                  # next biomass class x

}                                   # next time t

par(mfrow=c(1,1))

image(x=times,y=x_class[-1],z=bestpatch,ylab="weight",xlab="time",zlim=c(0,3),
      main="optimal patch",col=c("black","darkgrey","lightgrey","white"))
box()
legend("topleft",fill=c("black","darkgrey","lightgrey","white"),legend=c("dead","1","2","3"))
contour(x=1:t_max,y=x_class,z=f,add=TRUE)
legend("topright",legend="fitness",lty=1)


subtitle()

par(ask=opar$ask)
par(mar=opar$mar)
par(oma=opar$oma)
par(mfrow=opar$mfrow)

