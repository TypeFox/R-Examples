cosvol=function(area=60, zmax=1, zmin=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, inunit='deg2', ref){
  if(!missing(ref)){
    params=.getcos(ref)
    H0=params['H0']
    OmegaM=params['OmegaM']
    OmegaL=params['OmegaL']
  }
  if(inunit=='amin2'){area=area/3600}
  if(inunit=='asec2'){area=area/12960000}
  if(inunit=='rad2' | inunit=='sr'){area=area*(180/pi)^2}
  vols=cosdistCoVol(z=c(zmin,zmax), H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL)
  totalvol=as.numeric((vols[2]-vols[1])*area*pi/129600)
  volmeanz=0.75*(zmax^4-zmin^4)/(zmax^3-zmin^3)
  volmedz=((zmax^3+zmin^3)/2)^(1/3)
  return(c(voltot=totalvol, volmeanz=volmeanz, volmedz=volmedz))
}