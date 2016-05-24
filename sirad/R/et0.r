et0 <-
  function(Tmax,Tmin, vap_pres,sol_rad,tal,z,uz,meah=10,extraT=NA,days=NA,lat=NA) {
    
    if (!(is.numeric(extraT) & length(extraT)==length(Tmax)))  {
      extraT <- extrat(dayOfYear(days), radians(lat))$ExtraTerrestrialSolarRadiationDaily  
    }
    Tmean <- (Tmax+Tmin)/2
    
    et0 <- (0.408*deltaVP(Tmax,Tmin)*(rns(sol_rad,albedo=0.23)-rnl(Tmax,Tmin,sol_rad,vap_pres,extraT,tal))+psychC(Tmax,Tmin,z)*(900/(Tmean+273))*wind2(uz,meah)*(es(Tmax,Tmin)-vap_pres))/(deltaVP(Tmax,Tmin)+psychC(Tmax,Tmin,z)*(1+0.34*wind2(uz,meah)))
    et0
  }