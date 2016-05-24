#'@name k.read
#'@aliases 
#'k.read
#'k.read.soloviev
#'k.cole
#'k.macIntyre
#'k.crusius
#'k.vachon
#'k.heiskanen
#'@title Returns a timeseries of gas exchange velocity
#'@description 
#'Returns the gas exchange velocity based on the chosen model in units of m/day
#'@usage
#'k.cole(ts.data)
#'
#'k.crusius(ts.data, method='power')
#'
#'k.read(ts.data, wnd.z, Kd, atm.press, lat, lake.area)
#'
#'k.read.soloviev(ts.data, wnd.z, Kd, atm.press, lat, lake.area)
#'
#'k.macIntyre(ts.data, wnd.z, Kd, atm.press,params=c(1.2,0.4872,1.4784))
#'
#'k.vachon(ts.data, lake.area, params=c(2.51,1.48,0.39))
#'
#'k.heiskanen(ts.data, wnd.z, Kd, atm.press)
#'@param ts.data vector of datetime in POSIXct format
#'@param method Only for \link{k.crusius}. String of valid method . Either "linear", "bilinear", or "power"
#'@param wnd.z height of wind measurement, m
#'@param Kd Light attenuation coefficient (Units:m^-1)
#'@param atm.press atmospheric pressure in mb
#'@param lat Latitude, degrees north
#'@param lake.area Lake area, m^2
#'@param params Only for \link{k.vachon.base} and \link{k.macIntyre}. See details.
#'
#'@details Can change default parameters of MacIntyre and Vachon models. Default for Vachon is 
#'c(2.51,1.48,0.39). Default for MacIntyre is c(1.2,0.4872,1.4784). Heiskanen 2014 uses MacIntyre 
#'model with c(0.5,0.77,0.3) and z.aml constant at 0.15.
#'
#'@return Returns a data.frame with a datetime column and a k600 column. k600 is in units of meters per day (m/d).
#'@import rLakeAnalyzer
#'@useDynLib LakeMetabolizer
#'@keywords methods math
#'
#'@references
#'Cole, J., J. Nina, and F. Caraco. \emph{Atmospheric exchange of carbon dioxide 
#'in a low-wind oligotrophic lake measured by the addition of SF~ 6}. 
#'Limnology and Oceanography 43 (1998): 647-656.
#'
#'MacIntyre, Sally, Anders Jonsson, Mats Jansson, Jan Aberg, Damon E. Turney, 
#'and Scott D. Miller. \emph{Buoyancy flux, turbulence, and the gas transfer 
#'coefficient in a stratified lake}. Geophysical Research Letters 37, no. 24 (2010).
#'
#'Read, Jordan S., David P. Hamilton, Ankur R. Desai, Kevin C. Rose, Sally MacIntyre, 
#'John D. Lenters, Robyn L. Smyth et al. \emph{Lake-size dependency of wind shear and convection 
#'as controls on gas exchange}. Geophysical Research Letters 39, no. 9 (2012).
#'
#'Crusius, John, and Rik Wanninkhof. \emph{Gas transfer velocities measured at low 
#'wind speed over a lake}. Limnology and Oceanography 48, no. 3 (2003): 1010-1017.
#'
#'Dominic Vachon and Yves T. Prairie. \emph{The ecosystem size and shape dependence 
#'of gas transfer velocity versus wind speed relationships in lakes}.
#'Can. J. Fish. Aquat. Sci. 70 (2013): 1757-1764.
#'
#'Jouni J. Heiskanen, Ivan Mammarella, Sami Haapanala, Jukka Pumpanen, Timo Vesala, Sally MacIntyre
#'Anne Ojala.\emph{ Effects of cooling and internal wave motions on gas
#'transfer coefficients in a boreal lake}. Tellus B 66, no.22827 (2014)
#'
#'Alexander Soloviev, Mark Donelan, Hans Graber, Brian Haus, Peter Schlussel.
#'\emph{An approach to estimation of near-surface turbulence and CO2 transfer
#'velocity from remote sensing data}. Journal of Marine Systems 66, (2007): 182-194.
#'
#'@author
#'Hilary Dugan, Jake Zwart, Luke Winslow, R. Iestyn. Woolway, Jordan S. Read
#'@seealso 
#'\link{k.cole}
#'\link{k.crusius}
#'\link{k.macIntyre}
#'\link{k.vachon}
#'\link{k.heiskanen}
#'@examples 
#'data.path = system.file('extdata', package="LakeMetabolizer")
#'
#'tb.data = load.all.data('sparkling', data.path)
#'
#'ts.data = tb.data$data #pull out just the timeseries data
#'
#'#calculate U10 and add it back onto the original 
#'
#'u10 = wind.scale(ts.data)
#'ts.data = rmv.vars(ts.data, 'wnd', ignore.offset=TRUE) #drop old wind speed column
#'ts.data = merge(ts.data, u10)                          #merge new u10 into big dataset  
#'
#'
#'k600_cole = k.cole(ts.data)
#'
#'k600_crusius = k.crusius(ts.data)
#'
#'kd        = tb.data$metadata$averagekd
#'wnd.z      = 10   #because we converted to u10
#'atm.press  = 1018
#'lat       = tb.data$metadata$latitude
#'lake.area = tb.data$metadata$lakearea
#'
#'#for k.read and k.macIntyre, we need LW_net. 
#'#Calculate from the observations we have available. 
#'
#'lwnet = calc.lw.net(ts.data, lat, atm.press)
#'ts.data = merge(ts.data, lwnet)
#'\dontrun{
#'k600_read = k.read(ts.data, wnd.z=wnd.z, Kd=kd, atm.press=atm.press, 
#'	lat=lat, lake.area=lake.area)
#'
#'k600_soloviev = k.read.soloviev(ts.data, wnd.z=wnd.z, Kd=kd, 
#'	atm.press=atm.press, lat=lat, lake.area=lake.area)
#'
#'k600_macIntyre = k.macIntyre(ts.data, wnd.z=wnd.z, Kd=kd, atm.press=atm.press)
#'}
#'@export
k.read = function(ts.data, wnd.z, Kd, atm.press, lat, lake.area){
  
  S_B <- 5.67E-8 # Stefan-Boltzman constant (K is used)
  emiss <- 0.972 # emissivity;
  Kelvin = 273.15 #conversion from C to Kelvin
  
  data = ts.data
  # Get short wave radiation data 
  if(has.vars(data, 'sw')){ 
    sw <- get.vars(data, 'sw')
    
  } else if (has.vars(data, 'par')){
    tmp.par = get.vars(data, 'par')
    sw = par.to.sw(tmp.par)
  } else {  
    stop("Data must have PAR or SW column\n")
  }
  
  wtr <- get.vars(data,'wtr')
  Ts <- get.Ts(data)

  airT <- get.vars(data, 'airt')
  
  RH <- get.vars(data, 'rh')
  
  # Get long wave radiation data
  if(has.vars(data, 'lwnet')){ 
    lwnet <- get.vars(data,'lwnet')
    
  } else if(has.vars(data, 'lw')){
    lw_in <- get.vars(data, 'lw') # long wave in
    Tk <- Ts+Kelvin # water temperature in Kelvin
    LWo <- S_B*emiss*Tk^4 # long wave out
    lwnet <- lw_in[,2]-LWo
    
    lwnet = data.frame(datetime=lw_in$datetime, lwnet=lwnet)
    
  } else {  
    stop("no longwave radiation available")
  }
  

  wnd <- get.vars(data, 'wnd')
  
  m.d = ts.meta.depths(wtr)
  
  k600 = k.read.base(wnd.z, Kd, lat, lake.area, atm.press, data$datetime, Ts[,2], m.d$top, 
                airT[,2], wnd[,2], RH[,2], sw[,2], lwnet[,2])
  
  return(data.frame(datetime=data$datetime, k600=k600))
}

#'@name k.read.base
#'@aliases 
#'k.read.base
#'k.read.soloviev.base
#'k.cole.base
#'k.macIntyre.base
#'k.crusius.base
#'k.vachon.base
#'k.heiskanen.base
#'@title Returns a timeseries of gas exchange velocity
#'@description 
#'Returns the gas exchange velocity based on the chosen model in units of m/day
#'@usage
#'k.cole.base(wnd)
#'
#'k.crusius.base(wnd, method='power')
#'
#'k.read.base(wnd.z, Kd, lat, lake.area, atm.press, dateTime, Ts, z.aml, 
#'airT, wnd, RH, sw, lwnet)
#'
#'k.read.soloviev.base(wnd.z, Kd, lat, lake.area, atm.press, dateTime, Ts, z.aml, 
#'airT, wnd, RH, sw, lwnet)
#'
#'k.macIntyre.base(wnd.z, Kd, atm.press, dateTime, Ts, z.aml, airT, wnd, RH, sw, 
#'lwnet, params=c(1.2,0.4872,1.4784))
#'
#'k.vachon.base(wnd, lake.area, params=c(2.51,1.48,0.39))
#'
#'k.heiskanen.base(wnd.z, Kd, atm.press, dateTime, Ts, z.aml, airT, wnd, RH, sw, lwnet)
#'
#'@param wnd Numeric value of wind speed, (Units:m/s)
#'@param method Only for \link{k.crusius.base}. String of valid method . Either "constant", "bilinear", or "power"
#'@param wnd.z Height of wind measurement, (Units: m)
#'@param Kd Light attenuation coefficient (Units: m^-1)
#'@param lat Latitude, degrees north
#'@param lake.area Lake area, m^2
#'@param atm.press Atmospheric pressure, (Units: millibar)
#'@param dateTime datetime (Y-\%m-\%d \%H:\%M), (Format: \code{\link{POSIXct}})
#'@param Ts Numeric vector of surface water temperature, (Units(deg C)
#'@param z.aml Numeric vector of actively mixed layer depths. Must be the same length as the Ts parameter
#'@param airT Numeric value of air temperature, Units(deg C)
#'@param RH Numeric value of relative humidity, \%
#'@param sw Numeric value of short wave radiation, W m^-2
#'@param lwnet Numeric value net long wave radiation, W m^-2
#'@param params Optional parameter input, only for \link{k.vachon.base} and \link{k.macIntyre.base}. See details.
#'@details Can change default parameters of MacIntyre and Vachon models. Default for Vachon is 
#'c(2.51,1.48,0.39). Default for MacIntyre is c(1.2,0.4872,1.4784). Heiskanen et al. (2014) uses MacIntyre 
#'model with c(0.5,0.77,0.3) and z.aml constant at 0.15.
#'@return Numeric value of gas exchange velocity (k600) in units of m/day. Before use, 
#'should be converted to appropriate gas using \link{k600.2.kGAS}.
#'@keywords methods math
#'@references
#'Cole, J., J. Nina, and F. Caraco. \emph{Atmospheric exchange of carbon dioxide 
#'in a low-wind oligotrophic lake measured by the addition of SF~ 6}. 
#'Limnology and Oceanography 43 (1998): 647-656.
#'
#'MacIntyre, Sally, Anders Jonsson, Mats Jansson, Jan Aberg, Damon E. Turney, 
#'and Scott D. Miller. \emph{Buoyancy flux, turbulence, and the gas transfer 
#'coefficient in a stratified lake}. Geophysical Research Letters 37, no. 24 (2010).
#'
#'Read, Jordan S., David P. Hamilton, Ankur R. Desai, Kevin C. Rose, Sally MacIntyre, 
#'John D. Lenters, Robyn L. Smyth et al. \emph{Lake-size dependency of wind shear and convection 
#'as controls on gas exchange}. Geophysical Research Letters 39, no. 9 (2012).
#'
#'Crusius, John, and Rik Wanninkhof. \emph{Gas transfer velocities measured at low 
#'wind speed over a lake}. Limnology and Oceanography 48, no. 3 (2003): 1010-1017.
#'
#'Dominic Vachon and Yves T. Prairie. \emph{The ecosystem size and shape dependence 
#'of gas transfer velocity versus wind speed relationships in lakes}.
#'Can. J. Fish. Aquat. Sci. 70 (2013): 1757-1764.
#'
#'Jouni J. Heiskanen, Ivan Mammarella, Sami Haapanala, Jukka Pumpanen, Timo Vesala, Sally MacIntyre
#'Anne Ojala.\emph{ Effects of cooling and internal wave motions on gas
#'transfer coefficients in a boreal lake}. Tellus B 66, no.22827 (2014)
#'
#'Alexander Soloviev, Mark Donelan, Hans Graber, Brian Haus, Peter Schlussel.
#'\emph{An approach to estimation of near-surface turbulence and CO2 transfer
#'velocity from remote sensing data}. Journal of Marine Systems 66, (2007): 182-194.
#'
#'@author
#'R. Iestyn. Woolway, Hilary Dugan, Luke Winslow, Jordan S Read, GLEON fellows
#'@seealso 
#'\link{k.cole}
#'\link{k.read}
#'\link{k.crusius}
#'\link{k.macIntyre}
#'\link{k.vachon}
#'\link{k.heiskanen}
#'@examples 
#'wnd.z <- 2
#'Kd <- 2
#'lat <- 54
#'lake.area <- 5000 
#'atm.press <- 1013
#'dateTime <- as.POSIXct("2013-12-30 14:00")
#'Ts <- 16.5
#'z.aml <- 2.32
#'airT <- 20
#'wnd <- 6
#'RH <- 90
#'sw <- 800
#'lwnet <- -55
#'timeStep <- 30
#'
#'U10 <- wind.scale.base(wnd, wnd.z)
#'
#'k600_cole <- k.cole.base(U10)
#'
#'k600_crusius <- k.crusius.base(U10)
#'
#'k600_read <- k.read.base(wnd.z, Kd, lat, lake.area, atm.press, 
#'dateTime, Ts, z.aml, airT, wnd, RH, sw, lwnet)
#'
#'k600_soloviev <- k.read.soloviev.base(wnd.z, Kd, lat, lake.area, 
#'atm.press, dateTime, Ts, z.aml, airT, wnd, RH, sw, lwnet)
#'
#'k600_macInytre <- k.macIntyre.base(wnd.z, Kd, atm.press, 
#'dateTime, Ts, z.aml, airT, wnd, RH, sw, lwnet)
#'
#'@export
k.read.base <- function(wnd.z, Kd, lat, lake.area, atm.press, dateTime, Ts, z.aml, airT, wnd, RH, sw, lwnet){ 
  
  # define constants used in function
  Kelvin <- 273.15 # temp mod for deg K   
  emiss <- 0.972 # emissivity;
  S_B <- 5.67E-8 # Stefan-Boltzman constant (?K is used)
  vonK <- 0.41 # von Karman  constant
  dT <- 0.5   # change in temp for mixed layer depth
  C1 <- 114.278 # from Soloviev et al. 2007
  nu <- 0.29 # proportionality constant from Zappa et al. 2007, lower bounds
  KeCrit <- 0.18     # constant for wave age = 20 (Soloviev et al. 2007)
  albedo_SW <- 0.07
  swRat <- 0.46 # percentage of SW radiation that penetrates the water column
  g <- 9.81 # gravity
  C_w <- 4186 # J kg-1 ?C-1 (Lenters et al. 2005)
  mnWnd <- 0.2 # minimum wind speed
  
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
    wnd <- wnd/(1-e1/vonK*log(10/wnd.z))
  }
  rhoAir <- 1.2 #  air density
  tau <- C_D*wnd^2*rhoAir
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
  
  # calculate buoyancy flux and w*
  B1 <- H_star*tExp*g
  B2 <- rho_w*C_w
  Bflx <- B1/B2
  ltI <- Bflx>0
  B1 <- Bflx
  B1[ltI] <- 0
  divi <- 1/3
  w1 <- -B1*z.aml
  wSt <- w1^divi
  
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
  
  KeDe <- (kinV*g)
  KeNm <- uSt^3
  Ke <- KeNm/KeDe
  tau <- tau    # estimate of total tau (includes wave stress)
  euPw <- (1+Ke/KeCrit)  # tau_shear = tau/(1+Ke/Kecr)
  tau_t <- tau/euPw      # tau_shear, Soloviev
  uTanS <- tau_t/rho_w   
  uTanS <- uTanS^0.5
  
  # calculate viscous sublayer
  Sv <- C1*kinV/uTanS
  eu_N <- uTanS^3      # e_u(0) = (tau_t/rho)^1.5/(vonK*Sv)
  eu_D <- vonK*Sv       # denominator
  eu_0 <- eu_N/eu_D    # in m2/s3
  ew_0 <- -1.0*B1       # buoyancy flux, but only when outward
  e_0 <- ew_0+eu_0     # e(0) from Soloviev (w/o wave effects)
  K1 <- e_0*kinV       # in units of m4/s4, want cm4/hr4
  K2 <- ew_0*kinV      # convective component (m4/s4)
  K1 <- K1*100^4*3600^4 # now in cm4/hr4  (Total)
  K2 <- K2*100^4*3600^4 # now in cm4/hr4  (Convective)
  K600 <- nu*600^(-0.5)*K1^(1/4)   # in cm/hr (Total)
  
  k600 <- as.numeric(K600)
  k600 <- k600*24/100 #now in units in m d-1
  return(k600)
}
