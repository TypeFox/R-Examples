

##################################################################
######  OMEXDIA: C, N, O2 diagenesis: simplified version    ######
###### a more flexible version will be published in         ######
###### a separate package - name: reacTran (to be decided)  ######
##################################################################

#====================#
# Model equations    #
#====================#

OMEXDIAsteady <- function () 
{
parms<- c(
# organic matter dynamics  #
MeanFlux = 20000/12*100/365,  # nmol/cm2/d - Carbon deposition: 20gC/m2/yr
rFast    = 0.01         ,  #/day        - decay rate fast decay detritus
rSlow    = 0.00001      ,  #/day        - decay rate slow decay detritus
pFast    = 0.9          ,  #-           - fraction fast detritus in flux
w        = 0.1/1000/365 ,  # cm/d       - advection rate
NCrFdet  = 0.16         ,  # molN/molC  - NC ratio fast decay detritus
NCrSdet  = 0.13         ,  # molN/molC  - NC ratio slow decay detritus

# oxygen and DIN dynamics  #

# Nutrient bottom water conditions
bwO2            = 300   ,  #mmol/m3     Oxygen conc in bottom water
bwNO3           = 10    ,  #mmol/m3
bwNH3           = 1     ,  #mmol/m3
bwODU           = 0     ,  #mmol/m3

# Nutrient parameters- inhibition and half-saturation cts changed compared to
# omexdia: 
NH3Ads          = 1.3   ,  #-           Adsorption coeff ammonium
rnit            = 500.  ,  #was:20 /d   Max nitrification rate
ksO2nitri       = 1.    ,  #umolO2/m3   half-sat O2 in nitrification
rODUox          = 20.   ,  #/d          Max rate oxidation of ODU
ksO2oduox       = 1.    ,  #mmolO2/m3   half-sat O2 in oxidation of ODU
ksO2oxic        = 3.    ,  #mmolO2/m3   half-sat O2 in oxic mineralisation
ksNO3denit      = 30.   ,  #mmolNO3/m3  half-sat NO3 in denitrification
kinO2denit      = 1.e-6 ,  #was:1. mmolO2/m3  half-sat O2 inhib denitrification
kinNO3anox      = 1.e-6 ,  #was:1. mmolNO3/m3 half-sat NO3 inhib anoxic degr
kinO2anox       = 1.e-6 ,  #was:1. mmolO2/m3  half-sat O2 inhib anoxic min

# Diffusion coefficients, temp = 10dgC
#Temp            = 10                     ,   # temperature 
DispO2          = 0.955    +10*0.0386    ,  #cm2/d
DispNO3         = 0.844992 +10*0.0336    ,
DispNH3         = 0.84672  +10*0.0336    ,
DispODU         = 0.8424   +10*0.0242)

# The grid
N      <- 100
dx     <- rep(0.1,100)
dx.int <- rep(0.1,101)
pormid <- rep(0.9,100)
porint <- rep(0.9,101)
Db     <- rep(1/365,101)

# First the steady-state condition
OC   <- rep(10,6*N)
DIA  <- steady.band(y=OC,func="omexdiamod",initfunc="initomexdia",
                   initpar=c(parms,dx,dx.int,
                   pormid,porint,Db),nspec=6,
                   dllname="ecolMod",nout=8,positive=TRUE)
steady <-DIA$y

# state variables rearranged from ordering per slice -> per spec
ii       <- as.vector(t(matrix(ncol=6,1:(6*N))))   
steady[ii] <-steady
return(list(steady=steady,precis=attr(DIA,"precis"),Solved=attr(DIA,"steady")))
}
