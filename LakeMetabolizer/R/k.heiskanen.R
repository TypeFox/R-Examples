# ---Author: Hilary Dugan, 2014-11-06 --- 
# -- References ----
##Jouni J. Heiskanen, Ivan Mammarella, Sami Haapanala, Jukka Pumpanen, Timo Vesala, Sally MacIntyre
#Anne Ojala. 2010. Effects of cooling and internal wave motions on gas
#transfer coefficients in a boreal lake.
#Tellus B 2014, 66, 22827, http://dx.doi.org/10.3402/tellusb.v66.22827 

# INPUTS;
# wnd.z: Height of anemometer
# Kd: Diffuse attenuation coefficient
# atm.press: Atmospheric pressure in hPa or Mb. If unknown use '1013' for sea level 
# dateTime: date and time vector in POSIXct
# Ts: numeric vector of surface water temperature, Units(deg C)
# z.aml: Numeric vector of actively mixed layer depths. Must be the same length as Ts parameter
# airT: vector of air temperatures in C
# wnd: vector of wind speed in m/s
# Rh: vector of relative humidity in %
# sw: vector of shortwave radiation in W/m2
# lwnet: vector of net longwave radiation. If missing, run calc.lw.net function first

# OUTPUT: Numeric value of gas exchange velocity (k600) in units of m/day. 
#Before use, should be converted to appropriate gas using k600.2.kGAS.

#'@export
k.heiskanen <- function(ts.data, wnd.z, Kd, atm.press){
  
  S_B <- 5.67E-8 # Stefan-Boltzman constant (K is used)
  emiss <- 0.972 # emissivity;
  Kelvin = 273.15 #conversion from C to Kelvin
  
  # Get short wave radiation data 
  if(has.vars(ts.data, 'sw')){ 
    sw <- get.vars(ts.data, 'sw')
    
  } else if (has.vars(ts.data, 'par')){
    
    tmp.par = get.vars(ts.data, 'par')
    sw = par.to.sw(tmp.par)
  } else {  
    stop("Data must have PAR or SW column\n")
  }
  
  # Get water temperature data
  wtr <- get.vars(ts.data, 'wtr')
  Ts <- get.Ts(ts.data)
  
  airT <- get.vars(ts.data, 'airt')
  
  RH <- get.vars(ts.data, 'rh')
  
  # Get long wave radiation data
  if(has.vars(ts.data, 'lwnet')){ 
    lwnet <- get.vars(ts.data,'lwnet')
    
  } else if(has.vars(ts.data, 'lw')){
    lw_in <- get.vars(ts.data, 'lw') # long wave in
    Tk <- Ts+Kelvin # water temperature in Kelvin
    LWo <- S_B*emiss*Tk^4 # long wave out
    lwnet <- lw_in[,2]-LWo
    
    lwnet = data.frame(datetime=lw_in$datetime, lwnet=lwnet)
    
  } else {  
    stop("no longwave radiation available")
  }
  
  wnd <- get.vars(ts.data, 'wnd')
  m.d = ts.meta.depths(wtr)
  
  k600 = k.heiskanen.base(wnd.z, Kd, atm.press, ts.data$datetime, Ts[,2], m.d$top, 
                          airT[,2], wnd[,2], RH[,2], sw[,2], lwnet[,2])
  
  return(data.frame(datetime=ts.data$datetime, k600=k600))
  
}
#'@export
k.heiskanen.base <- function(wnd.z, Kd, atm.press, dateTime, Ts, z.aml, airT, wnd, RH, sw, lwnet){
  
  #Constants
  S_B <- 5.67E-8 # Stefan-Boltzman constant (K is used)
  emiss <- 0.972 # emissivity
  Kelvin = 273.15 #conversion from C to Kelvin
  albedo_SW <- 0.07
  vonK <- 0.41 #von Karman constant
  swRat <- 0.46 # percentage of SW radiation that penetrates the water column
  mnWnd <- 0.2 # minimum wind speed
  g <- 9.81 # gravity
  C_w <- 4186 # J kg-1 ?C-1 (Lenters et al. 2005)
  
  # impose limit on wind speed
  rpcI <- wnd < mnWnd
  wnd[rpcI] <- mnWnd
  
  
  # calculate sensible and latent heat fluxes
  mm <- calc.zeng(dateTime,Ts,airT,wnd,RH,atm.press,wnd.z)
  C_D <- mm$C_D # drag coefficient for momentum
  E <- mm$alh # latent heat flux
  H <- mm$ash # sensible heat flux
  
  # calculate total heat flux
  dUdt <- sw*0.93 - E - H + lwnet
  Qo <- sw*(1-albedo_SW)*swRat
  
  # calculate water density
  rho_w <- water.density(Ts)
  
  # calculate u*
  if (wnd.z != 10) {
    e1 <- sqrt(C_D)
    u10 <- wnd/(1-e1/vonK*log(10/wnd.z))
  }else{
    u10 = wnd
  }

  rhoAir <- 1.2 #  air density
  vonK <- 0.41 # von Karman  constant
  tau <- C_D*u10^2*rhoAir
  uSt <- sqrt(tau/rho_w)
  
  # calculate the effective heat flux
  q1 <- 2-2*exp(z.aml*-Kd)
  q2 <- z.aml*Kd
  q3 <- exp(z.aml*-Kd)
  H_star <- dUdt-Qo*(q1/q2-q3) # Kim 1976
  
  # calculate the thermal expansion coefficient 
  thermalExpFromTemp <- function(Ts) {    
    V <- 1       
    dT <- 0.001 
    T1 <- water.density(Ts)
    T2 <- water.density(Ts+dT)
    V2 <- T1/T2
    dv_dT <- (V2-V)/dT
    alpha <- dv_dT
    return (alpha)
  }
  tExp <- thermalExpFromTemp(Ts)
  
  B1 = H_star*tExp*g #Imberger 1985: Effective heat flux * thermal expansion of water * gravity
  B2 = rho_w*C_w # mean density of the water column * specific heat of water at constant pressure
  Bflx = B1/B2
  
  wstar = ifelse(Bflx < 0,(-Bflx*z.aml)^(1/3),0) #penetrative convective velocity Heiskanen 2014 (Imberger 1985)
  Hk   = sqrt((0.00015*u10)^2 + (0.07*wstar)^2) 
  Hk   = Hk*100*3600 # Heiskanen's K in cm/hr
  Hk600 = Hk*600^(-0.5)
  k600 <- Hk600*24/100 #units in m d-1
  return(k600)
}



