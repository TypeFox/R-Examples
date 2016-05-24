.getcos=function(ref){
  cosref = NULL
  data('cosref',envir = environment())
  if(ref %in% cosref[,'Ref']==FALSE){stop('Provided ref name is not allowed, must be one of 737 / 137 / Planck / WMAP9 / WMAP7 / WMAP5 / WMAP3 / WMAP1 / Millennium / GiggleZ. See ?cosref for details.')}
  out=as.numeric(cosref[cosref[,'Ref']==ref,])
  names(out)=colnames(cosref)
  return(out)
}

cosdist=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, age=FALSE, ref, error=FALSE){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaK=1-OmegaM-OmegaL
  Einv = function(z, OmegaM, OmegaL, OmegaK) {1/sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL)}
  if(age){Einvz = function(z, OmegaM, OmegaL, OmegaK){1/(sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL) * (1 + z))}}
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    HubDist = (299792.458/H0)
    temp = suppressWarnings(integrate(Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000L))
    CoDist = HubDist * temp$value
    if(error){
      if(z>0){
        RelError = abs(temp$abs.error/temp$value)
      }else{
        RelError = 0
      }
    }
    
      if(OmegaK==0){
        CoDistTran = CoDist
        CoVol = ((4/3) * pi * CoDist^3)/1e9
      }else{
        if(OmegaK>0){
          CoDistTran = HubDist*(1/sqrt(OmegaK))*sinh(sqrt(OmegaK)*CoDist/HubDist)
          CoVol = ((4*pi*HubDist^3/(2*OmegaK))*((CoDistTran/HubDist)*sqrt(1+OmegaK*(CoDistTran/HubDist)^2)-(1/sqrt(abs(OmegaK)))*asinh(sqrt(abs(OmegaK))*(CoDistTran/HubDist))))/1e9
        }
        if(OmegaK<0){
          CoDistTran = HubDist*(1/sqrt(abs(OmegaK)))*sin(sqrt(abs(OmegaK))*CoDist/HubDist)
          CoVol = ((4*pi*HubDist^3/(2*OmegaK))*((CoDistTran/HubDist)*sqrt(1+OmegaK*(CoDistTran/HubDist)^2)-(1/sqrt(abs(OmegaK)))*asin(sqrt(abs(OmegaK))*(CoDistTran/HubDist))))/1e9
        }
      }
      
      a=1/(1+z)
      LumDist = (1+z)*CoDistTran
      AngDist = CoDistTran/(1+z)
      if(z>=0){DistMod = 5*log10(LumDist)+25}else{DistMod=NA}
      AngSize = AngDist*(pi/(180*60*60))*1000
      
      if (age) {
        HT = (3.08568025e+19/(H0*31556926))/1e9
        UniAge = HT*integrate(Einvz, 0, Inf, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000L)$value
        zAge = HT*integrate(Einvz, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000L)$value
      }
      if(error){
        if (age) {
          return = c(z = z, a = a, CoDist = CoDist, LumDist = LumDist, AngDist = AngDist, CoDistTran=CoDistTran, DistMod = DistMod, AngSize = AngSize, CoVol = CoVol, HubTime = HT, UniAgeNow = UniAge, UniAgeAtz = UniAge - zAge, TravelTime = zAge, RelError=RelError)
        }
        else {
          return = c(z = z, a = a, CoDist = CoDist, LumDist = LumDist, AngDist = AngDist, CoDistTran=CoDistTran, DistMod = DistMod, AngSize = AngSize, CoVol = CoVol, RelError=RelError)
        }
      }else{
        if (age) {
          return = c(z = z, a = a, CoDist = CoDist, LumDist = LumDist, AngDist = AngDist, CoDistTran=CoDistTran, DistMod = DistMod, AngSize = AngSize, CoVol = CoVol, HubTime = HT, UniAgeNow = UniAge, UniAgeAtz = UniAge - zAge, TravelTime = zAge)
        }
        else {
          return = c(z = z, a = a, CoDist = CoDist, LumDist = LumDist, AngDist = AngDist, CoDistTran=CoDistTran, DistMod = DistMod, AngSize = AngSize, CoVol = CoVol)
        }
      }
    }
    return(as.data.frame(t(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))))
}

cosdistz=function(z = 1){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  return(z)
}

cosdista=function(z = 1){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  return(1/(1 + z))
}

cosdistCoDist=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaK=1-OmegaM-OmegaL
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    Einv = function(z, OmegaM, OmegaL, OmegaK) {1/sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL)}
    HubDist = (299792.458/H0)
    CoDist = HubDist * integrate(Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000L)$value
    return=CoDist
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosdistCoDistTran=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaK=1-OmegaM-OmegaL
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    Einv = function(z, OmegaM, OmegaL, OmegaK) {1/sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL)}
    HubDist = (299792.458/H0)
    CoDist = HubDist * integrate(Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000L)$value
    if(OmegaK==0){
      CoDistTran = CoDist
    }else{
      if(OmegaK>0){
        CoDistTran = HubDist*(1/sqrt(OmegaK))*sinh(sqrt(OmegaK)*CoDist/HubDist)
      }
      if(OmegaK<0){
        CoDistTran = HubDist*(1/sqrt(abs(OmegaK)))*sin(sqrt(abs(OmegaK))*CoDist/HubDist)
      }
    }
    return=CoDistTran
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosdistLumDist=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaK=1-OmegaM-OmegaL
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    Einv = function(z, OmegaM, OmegaL, OmegaK) {1/sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL)}
    HubDist = (299792.458/H0)
    CoDist = HubDist * integrate(Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000L)$value
    if(OmegaK==0){
      CoDistTran = CoDist
    }else{
      if(OmegaK>0){
        CoDistTran = HubDist*(1/sqrt(OmegaK))*sinh(sqrt(OmegaK)*CoDist/HubDist)
      }
      if(OmegaK<0){
        CoDistTran = HubDist*(1/sqrt(abs(OmegaK)))*sin(sqrt(abs(OmegaK))*CoDist/HubDist)
      }
    }
    LumDist = (1+z) * CoDistTran
    return=LumDist
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosdistAngDist=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaK=1-OmegaM-OmegaL
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    Einv = function(z, OmegaM, OmegaL, OmegaK) {1/sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL)}
    HubDist = (299792.458/H0)
    CoDist = HubDist * integrate(Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000L)$value
    if(OmegaK==0){
      CoDistTran = CoDist
    }else{
      if(OmegaK>0){
        CoDistTran = HubDist*(1/sqrt(OmegaK))*sinh(sqrt(OmegaK)*CoDist/HubDist)
      }
      if(OmegaK<0){
        CoDistTran = HubDist*(1/sqrt(abs(OmegaK)))*sin(sqrt(abs(OmegaK))*CoDist/HubDist)
      }
    }
    AngDist = CoDistTran / (1 + z)
    return=AngDist
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosdistDistMod=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaK=1-OmegaM-OmegaL
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    Einv = function(z, OmegaM, OmegaL, OmegaK) {1/sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL)}
    HubDist = (299792.458/H0)
    CoDist = HubDist * integrate(Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000L)$value
    if(OmegaK==0){
      CoDistTran = CoDist
    }else{
      if(OmegaK>0){
        CoDistTran = HubDist*(1/sqrt(OmegaK))*sinh(sqrt(OmegaK)*CoDist/HubDist)
      }
      if(OmegaK<0){
        CoDistTran = HubDist*(1/sqrt(abs(OmegaK)))*sin(sqrt(abs(OmegaK))*CoDist/HubDist)
      }
    }
    if(z>=0){DistMod = 5*log10(CoDistTran*(1+z))+25}else{DistMod=NA}
    return=DistMod
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosdistAngSize=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaK=1-OmegaM-OmegaL
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    Einv = function(z, OmegaM, OmegaL, OmegaK) {1/sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL)}
    HubDist = (299792.458/H0)
    CoDist = HubDist * integrate(Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000L)$value
    if(OmegaK==0){
      CoDistTran = CoDist
    }else{
      if(OmegaK>0){
        CoDistTran = HubDist*(1/sqrt(OmegaK))*sinh(sqrt(OmegaK)*CoDist/HubDist)
      }
      if(OmegaK<0){
        CoDistTran = HubDist*(1/sqrt(abs(OmegaK)))*sin(sqrt(abs(OmegaK))*CoDist/HubDist)
      }
    }
    AngSize = (CoDistTran / (1 + z)) * (pi/(180 * 60 * 60)) * 1000
    return=AngSize
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosdistCoVol=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaK=1-OmegaM-OmegaL
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    Einv = function(z, OmegaM, OmegaL, OmegaK) {1/sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL)}
    HubDist = (299792.458/H0)
    CoDist = HubDist * integrate(Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000L)$value
    if(OmegaK==0){
      CoDistTran = CoDist
      CoVol = ((4/3) * pi * CoDist^3)/1e9
    }else{
      if(OmegaK>0){
        CoDistTran = HubDist*(1/sqrt(OmegaK))*sinh(sqrt(OmegaK)*CoDist/HubDist)
        CoVol = ((4*pi*HubDist^3/(2*OmegaK))*((CoDistTran/HubDist)*sqrt(1+OmegaK*(CoDistTran/HubDist)^2)-(1/sqrt(abs(OmegaK)))*asinh(sqrt(abs(OmegaK))*(CoDistTran/HubDist))))/1e9
      }
      if(OmegaK<0){
        CoDistTran = HubDist*(1/sqrt(abs(OmegaK)))*sin(sqrt(abs(OmegaK))*CoDist/HubDist)
        CoVol = ((4*pi*HubDist^3/(2*OmegaK))*((CoDistTran/HubDist)*sqrt(1+OmegaK*(CoDistTran/HubDist)^2)-(1/sqrt(abs(OmegaK)))*asin(sqrt(abs(OmegaK))*(CoDistTran/HubDist))))/1e9
      }
    }
    return=CoVol
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosdistUniAgeNow=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaK=1-OmegaM-OmegaL
  Einvz = function(z, OmegaM, OmegaL, OmegaK){1/(sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL) * (1 + z))}
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    HT = (3.08568025e+19/(H0 * 31556926))/1e9
    UniAge = HT * integrate(Einvz, 0, Inf, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000L)$value
    return=UniAge
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosdistUniAgeAtz=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z>=0)){stop('All z must be >=0')}
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaK=1-OmegaM-OmegaL
  Einvz = function(z, OmegaM, OmegaL, OmegaK){1/(sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL) * (1 + z))}
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    HT = (3.08568025e+19/(H0 * 31556926))/1e9
    UniAge = HT * integrate(Einvz, 0, Inf, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000L)$value
    zAge = HT * integrate(Einvz, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000L)$value
    return=UniAge-zAge
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosdistTravelTime=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaK=1-OmegaM-OmegaL
  Einvz = function(z, OmegaM, OmegaL, OmegaK){1/(sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL) * (1 + z))}
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    HT = (3.08568025e+19/(H0 * 31556926))/1e9
    zAge = HT * integrate(Einvz, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000L)$value
    return=zAge
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosdistHubTime=function(H0 = 100){
 return((3.08568025e+19/(H0 * 31556926))/1e9)
}

cosdistRelError=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaK=1-OmegaM-OmegaL
  temp = function(z, OmegaM, OmegaL, OmegaK) {
    Einv = function(z, OmegaM, OmegaL, OmegaK) {1/sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL)}
    temp = integrate(Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000L)
      if(z>0){
        RelError = abs(temp$abs.error/temp$value)
      }else{
        RelError=0
      }
    return=RelError
  }
  return(Vectorize(temp)(z = z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}


