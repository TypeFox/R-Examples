# ---Author: Hilary Dugan, 2013-10-20 --- 
# Last update: 2014-02-01 

#'@export
k.crusius = function(ts.data, method='power'){
  
  if(!has.vars(ts.data, 'wnd')){
    stop('k.crusius requires a "wnd" (wind speed) column in the supplied data')
  }
  
  wind = get.vars(ts.data, 'wnd')
  
  k600 = k.crusius.base(wind[,2], method)
  
  return(data.frame(datetime=ts.data$datetime, k600=k600))
}


#'@export
k.crusius.base <- function(wnd, method='power'){

  U10 = wnd  #This function uses just the wind speed it is supplied. 
  method = tolower(method)
  if (method=="constant"){
    k600 <- ifelse(U10<3.7,1,5.14*U10-17.9)
  } else if (method=="bilinear") {
    k600 <- ifelse(U10<3.7,0.72*U10,4.33*U10-13.3)
  } else if (method=="power") {
    k600 <- 0.228*U10^2.2+0.168 # units in cm h-1
  } else {
    stop('method must be one of three options {power, constant, and bilinear}')
  }
  
  k600 <- k600*24/100 #units in m d-1
  return(k600)
}


# -- References 
# CRUSIUS, JOHN, AND RIK WANNINKHOF. 2003
# Gas transfer velocities measured at low wind speed over a lake.
# Limnology and Oceanography. 48(3): 1010:1017.
