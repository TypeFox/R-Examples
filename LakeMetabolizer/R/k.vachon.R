# ---Author: Hilary Dugan, 2014-11-06 --- 
# -- References ----
#Dominic Vachon and Yves T. Prairie (2013)
#The ecosystem size and shape dependence of gas transfer velocity 
#versus wind speed relationships in lakes
#Can. J. Fish. Aquat. Sci. 70: 1757-1764 dx.doi.org/10.1139/cjfas-2013-0241

# INPUTS;
# wnd: wind speed 10 m height (m/s). Use wind.scale function to scale wind speed if necessary. 
# lake.area: lake area (m2)
# parA,parB,parC = option to change parameters from  
# default values of parameters parA=2.51, parB=1.48, parC=0.39

# OUTPUT: Numeric value of gas exchange velocity (k600) in units of m/day. 
#Before use, should be converted to appropriate gas using k600.2.kGAS.

#'@export
k.vachon <- function(ts.data, lake.area, params=c(2.51,1.48,0.39)){
  if(!has.vars(ts.data, 'wnd')){
    stop('k.vachon requires a "wnd" column in the supplied data')
  }
  wnd <- get.vars(ts.data, 'wnd')
  k600 <- k.vachon.base(wnd[,2],lake.area,params)
  return(data.frame(datetime=ts.data$datetime, k600=k600))
}

#'@export
k.vachon.base <- function(wnd, lake.area, params=c(2.51,1.48,0.39)){
  U10 <- wnd  #This function uses just the wind speed it is supplied
  k600 <- params[1] + params[2]*U10 + params[3]*U10*log10(lake.area/1000000) # units in cm h-1
  k600 <- k600*24/100 #units in m d-1
  return(k600)
}

