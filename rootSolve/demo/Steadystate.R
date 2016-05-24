mf <- par("mfrow")
par (mfrow=c(1,1))
example(steady.band)
par (mfrow=mf)
example(steady.1D)


##################################################################
######  EXAMPLE:  NPZ riverine model                        ######
##################################################################
# Example from the book:
# Soetaert and Herman (2009).
# a practical guide to ecological modelling
# using R as a simulation platform
# Springer

#====================#
# Model equations    #
#====================#

NPZriver <- function (time=0,state,parms=NULL)
{
 N<- state[1       :Nb    ]      # nutrients
 P<- state[(Nb+1)  :(2*Nb)]      # phytoplankton
 Z<- state[(2*Nb+1):(3*Nb)]      # zooplankton

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
 dN    <- -diff(NFlux)/delx - Pprod + Graz*(1-eff) + Zmort
 dP    <- -diff(PFlux)/delx + Pprod - Graz 
 dZ    <- -diff(ZFlux)/delx         + Graz*eff     - Zmort 

return(list(c(dN=dN,dP=dP,dZ=dZ), # rate of changes
             Pprod=Pprod,         #Profile of primary production 
             Graz=Graz,           #Profile of zooplankton grazing 
             Zmort=Zmort,         #Profile of zooplankton mortality 
             Nefflux=NFlux[Nb],   # DIN efflux
             Pefflux=PFlux[Nb],   # Phytoplankton efflux
             Zefflux=ZFlux[Nb]))  # Zooplankton efflux
}

#====================#
# Model run          #
#====================#
# model parameters
riverlen <- 100               # total length river,         km
Nb       <- 100               # number boxes
delx     <- riverlen / Nb     # box length

RiverN   <- 100               # N concentration at river,   mmolN/m3
RiverP   <- 10                # P concentration at river,   mmolN/m3
RiverZ   <- 1                 # Z concentration at river,   mmolN/m3

flow     <- 1                 # river flow,                 km/day

mumax    <- 0.5               # maximal light-limited primary prod, /day
kN       <- 1                 # half-saturated N for pprod, mmolN/m3
gmax     <- 0.5               # max grazing rate,           /day
kP       <- 1                 # half-saturated P for graz , mmolN/m3
mrt      <- 0.05              # mortality rate,             /day
eff      <- 0.7               # growth efficiency,          -

# initial guess of N, P, Z
Conc     <- rep(1,3*Nb)  

# steady-state solution, v= 1 km/day
flow     <- 1                  # river flow,                 km/day
sol      <- steady.1D (Conc, func=NPZriver,nspec=3,maxiter=100,
                       parms=NULL,atol=1e-10,positive=TRUE) 
conc     <- sol$y

# steady-state solution, v=5 km d-1
flow     <- 5                  # river flow,                 km/day
Conc     <- rep(1,3*Nb)  
sol2     <- steady.1D (Conc, func=NPZriver,nspec=3,maxiter=100,
                       parms=NULL,atol=1e-10,positive=TRUE) 
conc     <- cbind(conc,sol2$y)

# steady-state solution, v=10 km d-1
flow     <- 10                  # river flow,                 km/day
Conc     <- rep(1,3*Nb)  
sol3  <- steady.1D (Conc, func=NPZriver,nspec=3,maxiter=100,
                    parms=NULL,atol=1e-10,positive=TRUE) 
conc  <- cbind(conc,sol3$y)

# rearranging output
N    <- conc[1       :Nb    ,]
P    <- conc[(Nb+1)  :(2*Nb),]
Z    <- conc[(2*Nb+1):(3*Nb),]

#====================#
# plotting 
#====================#
par(mfrow=c(2,2))
Dist <- seq(delx/2,riverlen,by=delx)
matplot(Dist,N,ylab="mmolN/m3",main="DIN",type="l",
        lwd=2,ylim=c(0,110),xlab="Distance, km")
matplot(Dist,P,ylab="mmolN/m3" ,main="Phytoplankton",type="l",
        lwd=2,ylim=c(0,110),xlab="Distance, km")
matplot(Dist,Z,ylab="mmolN/m3" ,main="Zooplankton",type="l",
        lwd=2,ylim=c(0,110),xlab="Distance, km")

plot(0,type="n",xlab="",ylab="",axes=FALSE)
legend("center",lty=c(1,1,2),lwd=c(2,1,1),col=1:3,
        title="river flowrate",
       c(expression(v==1~km~d^{-1}),expression(v==5~km~d^{-1}),
         expression(v==10~km~d^{-1}))) 

mtext(side=3,outer=TRUE,"NPZ model in a river",line=-1.5,cex=1.5)	


par(mfrow=mf)
