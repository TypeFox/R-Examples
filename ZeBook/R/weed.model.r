################################################################################
# Working with dynamic models for agriculture
# R script for practical work
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013-06-09
# Model described in the book, Appendix. Models used as illustrative examples: description and R code
################################ FUNCTIONS #####################################
#' @title The Weed model - calculate change for one year
#' @description TOCOMPLETE
#' @param d : weed density at seed emergence (plants/m2) - value for year
#' @param S : seed production per m2 - value for year
#' @param SSBa : surface seedbank after tillage (grains/m2) - value for year
#' @param DSBa : deep seedbank after tillage (grains/m2) - value for year
#' @param Soil : value for soil decision (1 or 0)
#' @param Crop : value for crop decision (1 or 0)
#' @param Herb : value for herbicide treatment decision (1 or 0)
#' @param param : vector of the 16 parameters
#' @return a vector with values of state variables for year+1
#' @export
weed.update=function(d, S, SSBa, DSBa, Soil, Crop, Herb, param) {
# 4 state variables
# S : seed production (seeds/m2)
# d : weed density at seed emergence (plants/m2)
# SSBa : surface seedbank after tillage (grains/m2)
# DSBa : deep seedbank after tillage (grains/m2)
# Soil, Crop, Herb : decisions
# param : 16 parameter values of the model, read from the param vector
# Intermediate variables

beta=param["beta.1"]*Soil+param["beta.0"]*(1-Soil)
chsi=param["chsi.1"]*Soil+param["chsi.0"]*(1-Soil)
Smax=param["Smax.1"]*Crop+param["Smax.0"]*(1-Crop)
alpha=Smax/160000

# Update State variable
# Evolution of weed grain stock (death, germination)
# in surface soil
SSBb1=(1-param["mu"])*(SSBa-d)+param["v"]*(1-param["phi"])*S
# in depth soil
DSBb1=(1-param["mu"])*DSBa

# Soil tillage effect
SSBa1=(1-beta)*SSBb1+chsi*DSBb1
DSBa1=(1-chsi)*DSBb1+beta*SSBb1

#d1=weed density at seed emergence (plants per m2) in year+1
d1=param["delta.new"]*param["v"]*(1-param["phi"])*(1-beta)*S + param["delta.old"]*(SSBa1-S*param["v"]*(1-param["phi"])*(1-beta))
#dm1=weed density at maturity depending on herbicide  (plants per m2) in year+1 (an auxillary variable). 
dm1=(1-param["mh"]*Herb)*(1-param["mc"])*d1
#S=seed production per m2 in year+1
S1=Smax*dm1/(1+alpha*dm1)
# Yield in year+1 (an auxillary variable).
Yield1=max(0,param["Ymax"]*(1-(param["rmax"]*dm1/(1+param["gamma"]*dm1))))

return(c(d1, S1, SSBa1, DSBa1, Yield1))
}
################################################################################
#' @title The Weed model - calculate daily values over designated time period
#' @description TOCOMPLETE
#' @param param : vector of the 16 parameters
#' @param weed.deci : decisions table for Soil, Crop et Herbicide 
#' @return data.frame with annual values of yield
#' @export
weed.model=function(param, weed.deci) {
  # Duration of the simulation according to weed.deci (years)
  duration=length(weed.deci$Soil)

  # Initialize variables
  # 4 states variables, as vectors initialized to NA
  # Weed density (plants/m2)
  d=rep(NA,duration+1)
  # Seed production (grains/m2)
  S=rep(NA,duration+1)
  # Surface seed bank after tillage (grains/m2)
  SSBa=rep(NA,duration+1)
  # Surface seed bank after tillage (grains/m2
  DSBa=rep(NA,duration+1)

  # Yield (ton/ha), an auxillary variable. 
  Yield=rep(NA,duration+1)

  # Initialize state variables when sowing on day "sdate"
  d[1]=400
  S[1]=68000
  SSBa[1]=3350
  DSBa[1]=280
  # compute Yield for initial year
  D0=(1-param["mh"]*weed.deci$Herb[1])*(1-param["mc"])*d[1]
  Yield[1]=max(0,param["Ymax"]*(1-(param["rmax"]*D0/(1+param["gamma"]*D0))))

  # Integration loop
  for (year in 1:(duration)) {
    # Update state variables
    Z=weed.update(d[year], S[year], SSBa[year], DSBa[year],weed.deci$Soil[year], weed.deci$Crop[year], weed.deci$Herb[year], param)
    d[year+1]=Z[1]
    S[year+1]=Z[2]
    SSBa[year+1]=Z[3]
    DSBa[year+1]=Z[4]
    Yield[year+1]=Z[5]
}
  # End simulation loop
#return(Yield)
return(data.frame(year=0:duration,d=d,S=S,SSBa=SSBa,DSBa=DSBa,Yield=Yield))
}

################################################################################
#' @title Wrapper function to run the Weed model multiple times (for multiple sets of inputs)
#' @param X : parameter matrix
#' @param weed.deci : decisions table for Soil, Crop et Herbicide 
#' @return matrix with Yield for year 3 for each parameter vector
#' @export
weed.simule=function(X, weed.deci) {
  # 4th row for year 3
  Y = apply(X,1,function(param) weed.model(param, weed.deci)[4,"Yield"])
  return(as.matrix(Y))
}

################################################################################
#' @title Define parameter values of the Weed model
#' @return matrix with parameter values (nominal, binf, bsup)
#' @export
weed.define.param = function()
{
# Nominal values of parameters
mu=0.84
v=0.6
phi=0.55
beta.1=0.95
beta.0=0.2
chsi.1=0.3
chsi.0=0.05
delta.new=0.15
delta.old=0.3
mh=0.98
mc=0
Smax.1=445
Smax.0=296
Ymax=8
rmax=0.002
gamma=0.005

param=data.frame("mu"=mu,"v"=v,"phi"=phi,"beta.1"=beta.1,"beta.0"=beta.0,
"chsi.1"=chsi.1,"chsi.0"=chsi.0,"delta.new"=delta.new,"delta.old"=delta.old,
"mh"=mh,"mc"=mc,"Smax.1"=Smax.1,"Smax.0"=Smax.0,"Ymax"=Ymax,"rmax"=rmax,
"gamma"=gamma, row.names = "nominal")


# definition of inferior and superior values.
param["binf",]=0.90*param["nominal",]
param["bsup",]=1.1*param["nominal",]
# for these 3 parameters, according to the equations, we need to specify bsup as :
param["bsup","mc"]=0.1
param["bsup","mh"]=min(1,param["bsup","mh"])
param["bsup","beta.1"]=min(1,param["bsup","beta.1"])
# nominal, binf, bsup
#row.names(param)=c("nominal","binf","bsup")
return(as.matrix(param))
}

# End file
