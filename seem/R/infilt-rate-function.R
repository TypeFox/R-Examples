# infiltration rate function
# infilt rate capacity based on deficit
# and limited by dry and sat infiltration
# infiltration capacity Dinkin & Nazimov (1995)

infilt.rate <- function(param, rain, v.soil){

  # infiltration rate parameters in mm/hr  
  Kd <- param[1]; Ks <- param[2]
  # soil saturation capacity in mm
  Z   <- param[3]
  # field capacity and  porosity
  Fc.coeff <- param[4]; porosity <- param[5]
  # soil sat capcity and field capacity 
  cap.soil <- Z * porosity
  Fc <- cap.soil * Fc.coeff

  # soil water deficit in mm
  # using deficit to determine infiltration cap
  deficit <- cap.soil - v.soil
  # base infilt is saturation infilt rate in mm/hr
  base <- Ks
  # slope or linear component mm/hr per mm
  linear <- (Kd-Ks)/cap.soil
  # infiltration rate in mm/hr
  infilt <- base + deficit*linear

  if (infilt > Kd) infilt <- Kd
  if (v.soil >= cap.soil) f<- Ks
  else   f<- infilt

  if(rain< f ) q <- rain else q <- f

  if(v.soil >= Fc) g <- Ks * (v.soil-Fc)/(cap.soil-Fc)
  else   g <- 0

 return(list(infilt.cap=f, infilt.rate=q, percol=g))
 # all in mm/hr		
}


