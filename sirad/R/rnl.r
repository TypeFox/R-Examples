rnl <-
  function(Tmax,Tmin,sol_rad,vap_pres,extraT,tal) {
    rnl <- 4.903e-09*(((Tmax+273.16)^4+(Tmin+273.16)^4)/2)*(0.34-0.14*sqrt(vap_pres))*(1.35*sol_rad/(extraT*tal)-0.35)
    rnl
  }