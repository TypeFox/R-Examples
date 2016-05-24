 # Functions return:
#
# tAdapt15251 - adaptive comfort temperature according to EN 15251
# tnAuliciems - neutral temperature based on Auliciems
# tnhumphreys - neutral temperature based on Humphreys
# mixR - mixing ratio of water in dry air [g/kg dry air]
# enth - enthalpy of air [kJ/kg]
# dewp - dewpoint of air [degree C]
# Humx - humidex of air [degree C]
# VP - current vapor pressure of air [pa]
# SVP - saturation vapor pressure of air [pa]
# rh - relative humidity [%] if only air temperature, mixIng ratio and barometric pressure are known

# v1.0 done by Marcel Schweiker, Michael Kleber, Sophia Mueller


# Function: Other models ################
###########################################
  
# Functions: adaptive comfort functions
calctAdapt15251 <- function(trm=20){
	data.frame(tAdapt15251 = 0.33*trm + 18.8)
	}
	
calctAdaptASHRAE <- function(tmmo){
	data.frame(tAdaptASHRAE = 0.33*tmmo + 17.8)
	}

calctnAuliciems <- function(ta, tmmo){
	data.frame(tnAuliciems = 9.22+0.48*ta+0.14*tmmo)
 } 

calctnHumphreysNV <- function(tmmo){ 
  data.frame(tnHumphreysNV = .534*tmmo + 11.9)
  }

calctnHumphreysAC <- function(tmmo){ 
  data.frame(tnHumphreysAC = 23.9+.295*(tmmo-22)*exp(-((tmmo-22)/(24*2 ^ .5)) ^ 2))
  }

# Functions: Other air humidity values ################
###########################################  

# mixIng ratio of water in dry air
calcMixR <- function(ta, rh, pb){
  1000*((exp((17.62*ta)/(243.12+ta))*611.2*(rh/100))/(461.51*(ta+273.15)))/(((pb/760*101325)-(exp((17.62*ta)/(243.12+ta))*611.2*(rh/100)))/(287.058*(ta+273.15)))
  }

# enthalpy of air
calcEnth <- function(ta, rh, pb){
  1.006*ta + (calcMixR(ta, rh, pb)/1000)*(1.86*ta + 2500)
  }
  
# dewpoint of air
calcDewp <- function(ta, rh){
  (rh/100) ^ (1/8.02)*(109.8+ta)-109.8
  }

# humidex of air
calcHumx <- function(ta, rh){
	ta+5/9*(6.11*exp(5417.753*(1/273.15-1/(calcDewp(ta, rh)+273.15)))-10)
  }

# saturation vapor pressure
calcSVP <- function(ta){
	6.1078*10^((7.5*ta)/(237.3+ta))*100
  }

# vapor pressure
calcVP <- function(ta, mr, pb){
	(mr/1000*(pb/760*101325)*462.51*(ta+273.15))/(mr/1000*(ta+273.15)*462.51+(ta+273.15)*287.058)
  }

# relative humidity from air temperature, mixIng ratio and barometric pressure
calcRH <- function(ta, mr, pb){
	calcVP(ta, mr, pb)/calcSVP(ta)*100
  }

