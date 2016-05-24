cosgrow=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, Sigma8=0.8, fSigma8=FALSE, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  OmegaK=1-OmegaM-OmegaL
  temp=function(z, H0, OmegaM, OmegaL, OmegaK){
    OmegaSum=OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL
    Hz=H0*sqrt(OmegaSum)
    OmegaMAtz=(OmegaM*(1+z)^3)/OmegaSum
    OmegaLAtz=OmegaL/OmegaSum
    OmegaKAtz=(OmegaK*(1+z)^2)/OmegaSum
    OmegaK=1-OmegaM-OmegaL
    Einva3=function(a, OmegaM, OmegaL, OmegaK){1/(a^3*(sqrt(OmegaM*a^(-3) + OmegaK*a^(-2) + OmegaL))^3)}
    Factor=(5*OmegaM/2)*(Hz/H0)*(1+z)*integrate(Einva3,0,1/(1+z),OmegaM=OmegaM,OmegaL=OmegaL,OmegaK=OmegaK,subdivisions=1000L)$value
    Factor0=(5*OmegaM/2)*integrate(Einva3,0,1,OmegaM=OmegaM,OmegaL=OmegaL,OmegaK=OmegaK,subdivisions=1000L)$value
    Sigma8Atz=Sigma8*(Factor/Factor0)/(1+z)
    if(fSigma8==FALSE){
      Rate=-1 - OmegaMAtz/2 + OmegaLAtz + (5*OmegaMAtz)/(2*Factor)
    }else{
      Rate=Sigma8Atz*(-1 - OmegaMAtz/2 + OmegaLAtz + (5*OmegaMAtz)/(2*Factor))
    }
    G=6.67384e-11 # m^3 kg^-1 s^-2
    km2m=1000
    Mpc2m=3.08567758e22
    Msol2kg=1.9891e30 # kg
    RhoCrit=(3*Hz)/(8*pi*G)*(km2m^2)*Mpc2m/Msol2kg #MsolperMpc3
  return=c(z=z, a=1/(1+z), H=Hz, OmegaM=OmegaMAtz, OmegaL=OmegaLAtz, OmegaK=OmegaKAtz, Factor=Factor, Rate=Rate, Sigma8=Sigma8Atz, RhoCrit=RhoCrit)
  }
  return(as.data.frame(t(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))))
}

cosgrowz=function(z = 1){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  return(z)
}

cosgrowa=function(z = 1){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  return(1/(1 + z))
}

cosgrowH=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaK=1-OmegaM-OmegaL
  return(H0*sqrt(OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL))
}

cosgrowOmegaM=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaK=1-OmegaM-OmegaL
  return((OmegaM*(1+z)^3)/(OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL))
}

cosgrowOmegaL=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaK=1-OmegaM-OmegaL
  return(OmegaL/(OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL))
}

cosgrowOmegaK=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaK=1-OmegaM-OmegaL
  return((OmegaK*(1+z)^2)/(OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL))
}

cosgrowFactor=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaK=1-OmegaM-OmegaL
  Einva3=function(a, OmegaM, OmegaL, OmegaK){1/(a^3*(sqrt(OmegaM*a^(-3) + OmegaK*a^(-2) + OmegaL))^3)}
  temp=function(z, OmegaM, OmegaL, OmegaK){
    growthfactor=(5*OmegaM/2)*cosgrowH(z,H0=1,OmegaM=OmegaM,OmegaL=OmegaL)*(1+z)*integrate(Einva3,0,1/(1+z),OmegaM=OmegaM,OmegaL=OmegaL,OmegaK=OmegaK,subdivisions=1000L)$value
    return=growthfactor
  }
  return(Vectorize(temp)(z = z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosgrowRate=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM, Sigma8=0.8, fSigma8=FALSE, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  OmegaMtemp=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  OmegaLtemp=cosgrowOmegaL(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  growthfacttemp=cosgrowFactor(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  if(fSigma8==FALSE){
    Sigma8temp=1
  }else{
    Sigma8temp=cosgrowSigma8(z=z, OmegaM=OmegaM, OmegaL=OmegaL,Sigma8=Sigma8)
  }
  return(Sigma8temp*(-1 - OmegaMtemp/2 + OmegaLtemp + (5*OmegaMtemp)/(2*growthfacttemp)))
}

cosgrowSigma8=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM, Sigma8=0.8, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  OmegaMtemp=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  OmegaLtemp=cosgrowOmegaL(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  growthfacttempAt0=cosgrowFactor(z=0, OmegaM=OmegaM, OmegaL=OmegaL)
  growthfacttempAtz=cosgrowFactor(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  return(Sigma8*(growthfacttempAtz/growthfacttempAt0)/(1+z))
}

cosgrowFactorApprox=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaMtemp=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  OmegaLtemp=cosgrowOmegaL(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  return((5*OmegaMtemp/2)/(OmegaMtemp^(4/7)-OmegaLtemp+(1+0.5*OmegaMtemp)*(1+OmegaLtemp/70)))
}

cosgrowRateApprox=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM, Sigma8=0.8, fSigma8=FALSE, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  if(fSigma8==FALSE){
    Sigma8temp=1
  }else{
    Sigma8temp=cosgrowSigma8(z=z, OmegaM=OmegaM, OmegaL=OmegaL,Sigma8=Sigma8)
  }
  OmegaMtemp=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  OmegaLtemp=cosgrowOmegaL(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  return(Sigma8temp*(OmegaMtemp^(4/7)+(1+OmegaMtemp/2)*(OmegaLtemp/70)))
}

cosgrowSigma8Approx=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM, Sigma8=0.8, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  OmegaMtemp=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  OmegaLtemp=cosgrowOmegaL(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  growthfacttempAt0=cosgrowFactorApprox(z=0, OmegaM=OmegaM, OmegaL=OmegaL)
  growthfacttempAtz=cosgrowFactorApprox(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  return(Sigma8*(growthfacttempAtz/growthfacttempAt0)/(1+z))
}

cosgrowRhoCrit=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaK=1-OmegaM-OmegaL
  G=6.67384e-11 # m^3 kg^-1 s^-2
  Hub2=cosgrowH(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL)^2 # (km/s / Mpc)^2
  km2m=1000
  Mpc2m=3.08567758e22
  Msol2kg=1.9891e30 # kg
  #rhocrit_kgperm3=(3*Hub2)/(8*pi*G) * ((km2m^2)/(Mpc2m^2)) #this is correct, should be ~9.2e-27 for H0=70 z=0
  #rhocrit_MsolperMpc3=rhocrit_kgperm3 * (Mpc2m^3)/Msol2kg #this is correct, should be ~2.8e11 for H0=100 z=0
  RhoCrit=(3*Hub2)/(8*pi*G) * (km2m^2)*Mpc2m/Msol2kg #compact form of the above, in units MsolperMpc3
  return(RhoCrit)
}
