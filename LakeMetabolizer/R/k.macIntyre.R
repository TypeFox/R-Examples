# ---Author: Hilary Dugan, 2013-10-20 --- 
# modified from K600 code by Richard Woolway and Jordan Read

# INPUTS;
# wnd.z: Height of anemometer
# Kd: Diffuse attenuation coefficient
# atm.press: Atmospheric pressure in hPa or Mb. If unknown use '1013' for sea level 
# dateTime: date and time vector
# wtr: dataframe of water temperatures 
# depth: vector of water temperature depths
# airT: vector of air temperatures in C
# Uz: vector of wind speed in m/s
# Rh: vector of relative humidity in %
# sw: vector of shortwave radiation in W/m2
# lw: vector of longwave radiation. If missing, run calc.lw.net function first
# par: vector of par data in umol m-2 s-1 

# OUTPUT: returns the gas exchange velocity for O2 in units of m/(timeStep*min) (i.e. 30 minute sampling 
#          interval will return kO2 in units of m/(1/48) - converts to fraction of day)

#'@export
k.macIntyre = function(ts.data, wnd.z, Kd, atm.press,params=c(1.2,0.4872,1.4784)){
  
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
  
  k600 = k.macIntyre.base(wnd.z, Kd, atm.press, ts.data$datetime, Ts[,2], m.d$top, 
                airT[,2], wnd[,2], RH[,2], sw[,2], lwnet[,2],params)
  
  return(data.frame(datetime=ts.data$datetime, k600=k600))
  
}
#'@export
k.macIntyre.base <- function(wnd.z, Kd, atm.press, dateTime, Ts, z.aml, airT, wnd, RH, sw, lwnet,
                             params=c(1.2,0.4872,1.4784)){
  
  #Constants
  S_B <- 5.67E-8 # Stefan-Boltzman constant (K is used)
  emiss <- 0.972 # emissivity;
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
  
  B1 = H_star*tExp*g
  B2 = rho_w*C_w
  Bflx = B1/B2
  
  
  # calculate kinematic viscosiy
  getKinematicVis <- function(Ts) {
    # from Mays 2005, Water Resources Engineering
    tempTable <- seq(0,100,by=5)
    # table in m2/s E-6
    visTable <- c(1.792,1.519,1.308,1.141,1.007,0.897,
                  0.804,0.727,0.661,0.605,0.556,0.513,0.477,0.444,
                  0.415,0.39,0.367,0.347,0.328,0.311,0.296)
    v <- data.frame(approx(tempTable,visTable,xout = Ts))[2]
    v <- v*1e-6
    return(v)
  }
  kinV <- getKinematicVis(Ts)
  KeNm = uSt^3
  
  #SmE   = 0.84*(-0.58*Bflx+1.76*KeNm/(vonK*z.aml))
  SmE = params[1]*-Bflx+params[2]*KeNm/(vonK*z.aml) #change to two coefficients
  SmE[SmE<0] = 0    # set negative to 0
  Sk   = SmE*kinV
  Sk   = Sk*100^4*3600^4 # Sally's K now in cm4/h4
  Sk600 = params[3]*600^(-0.5)*Sk^(1/4) # in cm/hr (Total) 
  
  k600 <- as.numeric(Sk600) # why is this not already numeric?
  k600 <- k600*24/100 #units in m d-1
  return(k600)
}

# -- References 
#MACINTRYE, sALLY, ANDERS JONSSON, MATS JANSSON, JAN ABERG, DAMON E. TURNEY AND SCOTT D. MILLER.
#2010. Buoyancy flux, turbulence, and the gas transfer coefficient in a stratified lake.
#Geophysical Research Letters. 37: L24604. doi:10.1029/2010GL044164 
