############################################################
# Example of model implementation in R and in compiled code
############################################################


require(rootSolve)

#==========================================================
# Parameters
#==========================================================

pars<- c(
D    = 1 ,    #/day, exchange rate with external world
Flux = 100,   # input flux of OM
r    = 0.1,   # consumption rate, /day
rox  = 1  ,   # HS oxidation rate, /day
ks   = 1  ,   # half saturation conc O2 limitation
ks2  = 1  ,   # half saturation conc SO4 limitation
BO2  = 100,   # boundary concentrations
BSO4 = 10000,
BHS  = 0)
y<-c(OM=1,O2=1,SO4=1,HS=1)


#==========================================================
# Model implemented in R
#==========================================================

model <- function (t,y,pars)
{
with (as.list(c(y,pars)),
{
  Min       = r*OM
  oxicmin   = Min*(O2/(O2+ks))
  anoxicmin = Min*(1-O2/(O2+ks))* SO4/(SO4+ks2)
  
  dOM  = Flux - oxicmin - anoxicmin
  dO2  = -oxicmin      -2*rox*HS*(O2/(O2+ks)) + D*(BO2-O2)
  dSO4 = -0.5*anoxicmin  +rox*HS*(O2/(O2+ks)) + D*(BSO4-SO4)
  dHS  =  0.5*anoxicmin  -rox*HS*(O2/(O2+ks)) + D*(BHS-HS)
  list(c(dOM,dO2,dSO4,dHS),SumS=SO4+HS)
})

}

# Another way to solve for steady-state:
# runs dynamically till state variables stop changing
#require(deSolve)
#root <- function(t,y,pars)
#{
# dy <- unlist(model(t,y,pars))
# sum(abs(dy))-1e-10
#}
#times<-c(0,1e16)
#print(system.time(
#for (i in 1:10)
#out<-lsodar(y=y,fun=model,times=times,parms=pars,rootfun=root)
#)/10 )


#-------------------------------------------------------------------------
# Steady-state solution, based on R-code
#-------------------------------------------------------------------------

print(system.time(
for (i in 1:100)
  ST <- steady(y = y, fun = model, parms = pars, pos = TRUE)
)/100)


#==========================================================
# Model implemented in Fortran
#==========================================================
# First compile anoxmod.f if not yet done.
# Make sure that the working directory is the directory with
# file anoxmod.f
# compiling is best done within R:
# system("R CMD SHLIB anoxmod.f") 
dyn.load("anoxmod.dll")
print(system.time(
for (i in 1:100)
  ST2 <- steady(y = y, fun = "anoxmod", parms = pars, dllname = "anoxmod",
              initfunc = "initanox", pos = TRUE, nout = 1)
)/100)
dyn.unload("anoxmod.dll")


#==========================================================
# Model implemented in C
#==========================================================
# First compile anoxmodc.c if not yet done.
# Make sure that the working directory is the directory with
# file anoxmodc.c
# compiling is best done within R:
# system("R CMD SHLIB anoxmodc.c") 
dyn.load("anoxmodc.dll")
print(system.time(
for (i in 1:100)
  ST3 <- steady(y = y, fun = "anoxmod", parms = pars, dllname = "anoxmodc",
              initfunc = "initanox", pos = TRUE, nout = 1)
)/100)
dyn.unload("anoxmodc.dll")


