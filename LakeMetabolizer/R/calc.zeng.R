
#'@title Estimate sensible and latent heat fluxes
#'@description 
#'Returns the sensible and latent heat fluxed based on Zeng et al, 1998'
#'@usage
#'calc.zeng(dateTime,Ts,airT,Uz,RH,atm.press,wnd.z,airT.z,RH.z)
#'
#'@param dateTime vector of datetime in POSIXct format
#'@param Ts numeric value of surface water temperature, degC
#'@param airT numeric value of air temperature, degC
#'@param Uz numeric value of wind speed, m/s
#'@param RH numeric value of relative humidity, \%
#'@param atm.press atmospheric pressure in mb
#'@param wnd.z height of wind measurement, m
#'@param airT.z height of air temperature measurement, m (optional)
#'@param RH.z height of relative humidity measurement, m (optional)
#'@return A data.frame including sensible and latent heat flux estimates, and other variables used in calculating these fluxes.
#'@keywords methods math
#'@importFrom stats complete.cases 
#'@references
#'Zeng, X., M. Zhao., and Dickinson, R.E. 1998. \emph{Intercomparison of bulk aerodynamic algorithms 
#'for the computation of sea surface fluxes using TOGA COARE and TAO data}. Journal of Climate 11: 2628-2644.
#'@author
#'R. Iestyn. Woolway
#'@seealso \link{k.read}
#'@examples 
#'dateTime <- as.POSIXct("2013-12-30 23:00")
#'Ts <- 22.51
#'airT <- 20
#'Uz <- 3  
#'RH <- 90
#'atm.press <- 1013
#'wnd.z <- 2
#'calc.zeng(dateTime,Ts,airT,Uz,RH,atm.press,wnd.z)
#'@export
calc.zeng <- function(dateTime,Ts,airT,Uz,RH,atm.press,wnd.z,airT.z,RH.z){
  
  psi <- function(k,zeta){
    chik <- (1 - 16*zeta)^0.25
    if (k == 1){
      psi <- 2*log((1 + chik)*0.5) + log((1+chik*chik)*0.5) - 2*atan(chik) + (pi/2)
    } else{
      psi <- 2*log((1 + chik*chik)*0.5)
    }
  }
  
  # generate data.frame of all variables (to ensure consistent times)
  dat <- data.frame(dateTime = dateTime,
                    Ts = Ts,
                    airT = airT,
                    Uz = Uz,
                    RH = RH)
  
  # remove duplicated time stamps
  dat$dateTime <- as.POSIXct(strptime(dat$dateTime,"%Y-%m-%d %H:%M")) # ensure times are POSIXct
  dat <- subset(dat,!duplicated(dat$dateTime)) #remove duplicate time stamps
  
  # store original dates - used for final data frame
  original_dates <- data.frame(dateTime = dat$dateTime)  
  
  # remove missing data
  dat <-  dat[complete.cases(dat),]
  
  # re-define variables
  Ts <- dat$Ts
  airT <- dat$airT
  Uz <- dat$Uz
  RH <- dat$RH
  
  # if temperature and humidity height are missing, assume same as wind
  if (missing(airT.z)){airT.z <- wnd.z}
  if (missing(RH.z)){RH.z <- wnd.z}

  # define constants
  const_vonKarman <- 0.41 # von Karman constant
  const_gas <- 287.1 # gas constant for dry air J kg-1 K-1
  const_SpecificHeatAir <- 1005 # Specific heat capacity of air, J kg-1 K-1
  const_Charnock <- 0.013 # charnock constant
  const_Gravity <- 9.81 # gravitational acceleration, m/s2

  # ensure wind below 0.2 are changed to 0.1
  thres <- 0.2
  Uz[Uz < thres] <- thres/2 # m/s

  # calculate humidity values
  e_s <- 6.11*exp(17.27*airT/(237.3 + airT)) # saturated vapour pressure at airT, mb
  e_a <- RH*e_s/100 # vapour pressure, mb
  q_z <- 0.622*e_a/atm.press # specific humidity, kg kg-1 
  e_sat <- 6.11*exp(17.27*Ts/(237.3 + Ts)) # saturated vapour pressure at Ts, mb
  q_s <- 0.622*e_sat/atm.press # humidity at saturation, kg kg-1

  # calculate gas constant for moist air, J kg-1 K-1
  R_a <- 287*(1 + 0.608*q_z)

  # -- calculate latent heat of vaporization, J kg-1 
  xlv <- 2.501e6-2370*Ts    

  # calculate air density, kg/m3
  rho_a <- 100*atm.press/(R_a*(airT + 273.16))

  # -- kinematic viscosity, m2 s-1
  KinV <- (1/rho_a)*(4.94e-8*airT + 1.7184e-5)   

  # calculate virtual air temperature, K
  t_virt <- (airT + 273.16)*(1 + 0.61*q_z)

  # estimate initial values of u* and zo
  ustar <- Uz*sqrt(0.00104+0.0015/(1+exp((-Uz+12.5)/1.56)))
  zo <- (const_Charnock*ustar^2/const_Gravity) + (0.11*KinV/ustar)
  for (i in 1:20){
    ustar <- const_vonKarman*Uz/(log(wnd.z/zo))
    zo <- (const_Charnock*ustar^2/const_Gravity) + (0.11*KinV/ustar)    
  }
  zo <- Re(zo)
  
  # calculate neutral transfer coefficients
  C_DN <- (ustar^2)/(Uz^2)
  re <- ustar*zo/KinV
  zot <- zo*exp(-2.67*(re)^(0.25) + 2.57)
  zot <- Re(zot)
  zoq <- Re(zot)
  C_HN <- const_vonKarman*sqrt(C_DN)/(log(wnd.z/zot)) 
  C_EN <- Re(C_HN)

  # calculate neutral transfer coefficients at 10 m
  C_D10N <- (const_vonKarman/log(10/zo))*(const_vonKarman/log(10/zo)) 
  C_E10N <- (const_vonKarman*const_vonKarman)/(log(10/zo)*log(10/zoq))
  C_H10N <- C_E10N    

  # calculate neutral latent and sensible heat fluxes, W/m2
  alhN <- rho_a*xlv*C_EN*Uz*(q_s-q_z)
  ashN <- rho_a*const_SpecificHeatAir*C_HN*Uz*(Ts-airT) 

  # calculate initial monin obukhov length scale, m
  obu <- (-rho_a*t_virt*(ustar*ustar*ustar))/(const_vonKarman*const_Gravity*(ashN/const_SpecificHeatAir + 0.61*(airT + 273.16)*alhN/xlv)) 

  # iteration to compute corrections for atmospheric stability    
  zeta_thres <- 15
  
  # pre-define arrays
  tstar <- rep(0,length(Uz))
  qstar <- rep(0,length(Uz))
  wc <- rep(0,length(Uz))
  
  for (i in 1:20){ 
    
    # calulate roughness lengths
    zo <- (0.013*((ustar^2)/const_Gravity)) + (0.11*(KinV/ustar))                       
    re <- (ustar*zo)/KinV     
    xq <- 2.67*re^0.25 - 2.57
    xq[xq < 0] <- 0
    zoq <- zo/exp(xq)        
    zot <- zoq     

    # define zeta
    zetam <- -1.574
    zetat <- -0.465
  
    # calculate ustar
    zeta <- wnd.z/obu
    
    zeta[zeta < -zeta_thres] <- -zeta_thres
    zeta[zeta > zeta_thres] <- zeta_thres
    
    idx <- zeta < zetam & !is.na(zeta) # very unstable conditions      
    ustar[idx] <- (Uz[idx]*const_vonKarman)/((log((zetam*obu[idx])/zo[idx]) - psi(1,zetam)) + 1.14*(((-zeta[idx])^0.333) - ((-zetam)^0.333))) 
    idx <- zeta < 0 & zeta >= zetam & !is.na(zeta) # unstable conditions
    ustar[idx] = (Uz[idx]*const_vonKarman)/(log(wnd.z/zo[idx]) - psi(1,zeta[idx]))
    idx <- zeta > 0 & zeta <= 1 & !is.na(zeta) # stable conditions
    ustar[idx] <- (Uz[idx]*const_vonKarman)/(log(wnd.z/zo[idx]) + 5*zeta[idx])
    idx <- zeta > 1 & !is.na(zeta) # very stable conditions
    ustar[idx] <- (Uz[idx]*const_vonKarman)/((log(obu[idx]/zo[idx])+ 5) + (5*log(zeta[idx]) + zeta[idx] - 1))
    
    # calculate tstar
    zeta <- airT.z/obu
    zeta[zeta < -zeta_thres] <- -zeta_thres
    zeta[zeta > zeta_thres] <- zeta_thres

    idx <- zeta < zetat & !is.na(zeta) # very unstable conditions
    tstar[idx] <- (const_vonKarman*(airT[idx] - Ts[idx]))/((log((zetat*obu[idx])/zot[idx]) - psi(2,zetat)) +  0.8*((-zetat)^-0.333 - ((-zeta[idx]))^-0.333))
    idx <- zeta >= zetat & zeta < 0 & !is.na(zeta) # unstable conditions
    tstar[idx] <- (const_vonKarman*(airT[idx] - Ts[idx]))/(log(airT.z/zot[idx]) - psi(2,zeta[idx]))
    idx <- zeta > 0 & zeta <= 1 & !is.na(zeta) # stable conditions
    tstar[idx] <- (const_vonKarman*(airT[idx] - Ts[idx]))/(log(airT.z/zot[idx]) + 5*zeta[idx])
    idx <- zeta > 1 & !is.na(zeta) # very stable conditions
    tstar[idx] <- (const_vonKarman*(airT[idx] - Ts[idx]))/((log(obu[idx]/zot[idx]) + 5) + (5*log(zeta[idx]) + zeta[idx] - 1))

    # calculate qstar
    zeta <- RH.z/obu
    zeta[zeta < -zeta_thres] <- -zeta_thres
    zeta[zeta > zeta_thres] <- zeta_thres

    idx <- zeta < zetat & !is.na(zeta) # very unstable conditions
    qstar[idx] <- (const_vonKarman*(q_z[idx] - q_s[idx]))/((log((zetat*obu[idx])/zoq[idx]) - psi(2,zetat)) + 0.8*((-zetat)^-0.333 - ((-zeta[idx]))^-0.333))
    idx <- zeta >= zetat & zeta < 0 & !is.na(zeta) # unstable conditions
    qstar[idx] <- (const_vonKarman*(q_z[idx] - q_s[idx]))/(log(RH.z/zoq[idx]) - psi(2,zeta[idx]))
    idx <- zeta > 0 & zeta <= 1 & !is.na(zeta) # stable conditions
    qstar[idx] <- (const_vonKarman*(q_z[idx] - q_s[idx]))/(log(RH.z/zoq[idx]) + 5*zeta[idx])
    idx <- zeta > 1 & !is.na(zeta) # very stable conditions
    qstar[idx] <- (const_vonKarman*(q_z[idx] - q_s[idx]))/((log(obu[idx]/zoq[idx]) + 5) + (5*log(zeta[idx]) + zeta[idx] - 1))
    
    # calculate zeta at 10 m
    zeta <- 10/obu
    zeta[zeta < -zeta_thres] <- -zeta_thres
    zeta[zeta > zeta_thres] <- zeta_thres

    # calculate transfer coefficients corrected for atmospheric stability
    C_D <- (ustar*ustar)/(Uz*Uz)

    # calculate tau and sensible and latent heat fluxes
    tau <- C_D*rho_a*ustar*ustar
    ash <- -rho_a*const_SpecificHeatAir*ustar*tstar
    alh <- -rho_a*xlv*ustar*qstar

    # calculate new monin obukhov length
    obu <- (-rho_a*t_virt*(ustar*ustar*ustar))/(const_Gravity*const_vonKarman*((ash/const_SpecificHeatAir) + (0.61*(airT + 273.16)*alh/xlv)))

    # alter zeta in stable cases
    zeta <- wnd.z/obu
    idx <- zeta >= 1
    Uz[idx] <- max(Uz[idx],0.1)

    # avoid singularity at um = 0 for unstable conditions     
    idx <- zeta < 0 & !is.na(zeta)
    th <- (airT + 273.16)*(1000/atm.press)^(287.1/1004.67) # potential temperature  
    thvstar <- tstar*(1 + 0.61*q_z/1000) + 0.61*th*qstar # temperature scaling parameter
    thv <- th*(1 + 0.61*q_z/1000) #virtual potential temperature    
  }

  # take real values to remove any complex values that arise from missing data or NA.
  C_D <- Re(C_D)
  
  # store results in data.frame and merge with original dateTime
  mm <- data.frame(dateTime = dat$dateTime,
                   C_D = C_D,alh = alh,ash = ash)
  mm <- merge(mm,original_dates,all = TRUE)
  return(mm)
  
}
