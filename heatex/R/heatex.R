# heatex.R
# created: 14 February 2013
# last revised: 18 February 2013


heatex<-function(df){

# define constants

ArAd<-0.70 # ratio of body area exposed to radiation versus the total body surface area. 0.70 for seated posture and 0.73 for standing
l<-2430  # latent heat of evaporation of sweat
m.coef<- 0.0001694 # permeance coefficient of the skin
emit<-0.97 # emittance from the outer surface of a clothed body 
sigma<-0.0000000567  # Stefan-Boltzmann constant
feff<-0.71 # effective radiation area of a clothed body

# Calculations

#--------------------------------
# Environmental variables
#---------------------------------

# convert deg C to Kelvin
tdb_K<-df$tdb + 273.16

# mean radiant temperature
tr<-((1+(0.222 * (df$va ^0.5))) * (df$tg - df$tdb)) + df$tdb

# convective heat transfer, hc
hc<-8.3 * (df$va * 0.6)

# radiative heat transfer coefficient, hr
hr<-4 * emit * sigma * ArAd * ((273.2 + ((df$tcl + tr)/2)) ^ 3)

# combined heat transfer coefficient, h
h<- hc + hr

# evaporative heat transfer coefficient, he
he<-16.5 * hc
  
# convert Pa from mmHg to kpa
pa_kpa<-df$pa * 0.1333

#---------------------------
# Physiological Variables
#---------------------------

# convert deg C to kelvin
tskf_K<-df$tskf + 273.16

# body surface area, AD
ad<-0.00718 * df$wt^0.425 * df$ht ^ 0.725

# mean body temperature, Tb

tbi<-(0.33 * df$tski + 0.67 * df$tci)
tbf<-(0.33 * df$tskf + 0.67 * df$tcf)

# saturated water vapor pressure at the skin surface, Ps

ps<-1.92 * df$tskf -25.3
ps_kpa<-ps * 0.1333 # convert mmHg to kPa
   
#----------------------------
# Clothing variables
#----------------------------

# clothing area factor, fcl
fcl<- 1 + (0.31 * (df$icl/0.155))

# effective clothing insulation, Icle
icle<-df$icl - ((fcl-1)/(0.155 * fcl * h))

# permeation efficiency factor of clothing, fpcl
fpcl<- 1/(1+(0.344 * hc * icle))

# intrinsic thermal resistance of clothing, Rc
rc<-(df$tskf - df$tdb)/hc

# intrinsic evaporative resistance of clothing, Re
re<-(ps_kpa - pa_kpa)/he

#-----------------------------------
# Partitional Calorimetry Equations
#-----------------------------------

# energy equivalent of oxygen, EE
ee<-(0.23 * df$rer + 0.77) * 21166

# metabolic free energy production, M
m<-(((ee * df$vo2 * df$time)/(df$time * 60))/ad)

# mechanical efficiency, n
n<-df$workrate/(m * ad)

# internal heat production, Hi
hi<-(m * (1-n))

# body heat storage, S
s<-((3474 * df$bmi * (tbf - tbi))/(df$time * 60))/ad

# heat transfer via conduction, K
# convert deg C to kelvin
tcl_K<-df$tcl + 273.16
k<-ad * ((tskf_K - tcl_K)/rc)

# heat transfer via radiation, R. For radiation from clothing surface, replace tskf with tcl.
r<-emit * sigma * fcl * feff *(df$tskf^4 - tr^4)

# heat transfer via convection, C. If convection from a clothed surface, change tskf_K to tcl_K.
conv<-(ad * fcl * hc * (tskf_K - tdb_K))/ ad

# required evaporative heat loss, Ereq
ereq<-hi - k - r- conv - s

# maximal evaporative capacity of the environment, Emax
emax<-fpcl * he * (ps_kpa - pa_kpa)

# skin wettedness, w
w<- ereq/ emax

# Evaporative heat transfer via skin diffusion, Ed
ed<-(l * m.coef * (ps - df$pa))

# Heat transfer by sweat evaporation from the skin surface, Esw
esw<-((((df$bmi*1000) - ((df$bmf*1000) + df$sweat - df$fluidfood - df$urinefaeces))-((0.019 * df$vo2 * (44 - df$pa)) * df$time))*2430)/((df$time * 60) * ad)

# Heat transfer via evaporation from the skin surface, Esk
esk<- ed + esw

# Return output
results<-data.frame(tr,hc,hr,h,he,pa_kpa,fcl,icle,fpcl,rc,re,ad,tbi,tbf,ps,ps_kpa,m,n,hi,s,k,r,conv,ereq,emax,w,ed,esw,esk)
    
} # End heatex function


