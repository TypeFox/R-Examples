# vector c attached to each node in the modwpt transform (Percival and Walden, page 215)
# vector c will be used to align LA wavelet coefficients
getC=function(j,n)
{
  vectorC=c();
  if (j==1)
  {
    if (n==0) vectorC=0
    else vectorC=1
    return(vectorC)
    
  }
  else
  {
    if (((n%%4)==0)||((n%%4)==3))
      return(c(getC(j-1,floor(n/2)),0))
    else 
      return(c(getC(j-1,floor(n/2)),1))
  }
  
}

# Auxiliary function used to calculate time shifts for wavelets(Percival and walden,229-230)
# This function will return c(Sjn0,Sjn1) 
getSc=function(j,n)
{
  vectorC=getC(j,n)
  Sjn1=(sum(vectorC*2^(0:(j-1))))
  Sjn0=2^j-1-Sjn1
  return (c(Sjn0,Sjn1))
  
  
}
# Auxiliary function used to calculate time shifts for LA wavelets(Percival and walden,229-230)
getLj=function(j,wf)
{
  #length of the wave filter
  L=length(wave.filter(wf)$lpf)
  Lj=(2^j-1)*(L-1)+1
  return(Lj)
}

#getVjnLA returns time shifts for la8,la16 and la20 wavelets in modwpt
getVjnLA=function(j,n,wf)
{
  aux=getSc(j,n);
  Sjn1=aux[2];
  return(abs(-getLj(j,wf)/2+(2^(j-1)-Sjn1)))
}

#performs left circular shift
circularShift=function(x,n)
{
  l=length(x)
  if ((n==0)||(n==l)) return(x)
  if ((n<0)||(n>l)) {n=n%%l}
  return(c(x[(n+1):l],x[1:n]))
}
# aligns <<x>> wavelet coefficients from modwpt using <<wf>> wavelet. <<x>> are the
# (j.n) node wavelet coefficients. <<wf>> should be "l8","l16" or "l20". If not, not 
# alignment is perform
align=function(x,j,n,wf)
{  
  #return(x)
  if ((wf!="la8")&&(wf!="la16")&&(wf!="la20")) 
  {    
    #general method
    return(circularShift(x,CenterOfEnergyAdvances(j,n,wf)))
  }
  else
  { #specific method for LA wavelets
    return(circularShift(x,getVjnLA(j,n,wf)))
  }
}


#shift applicable to (almost) every filter using the notion of the 'center of
#energy' by Hess-Nielsen and Wickerhauser.
# Walden page 231

# gets the appropiate advance for the (j.n) node of the <<wf>> wavelet. 
CenterOfEnergyAdvances=function(j,n,wf)
{
  #variables for the computation
  S=getSc(j,n)
  E=CenterEnergies(wf)
  advance=round(  S[1]*E[1]+S[2]*E[2]  ) 
  
  
  return (advance)
}

# gets the center of energy for both wf filters
CenterEnergies=function(wf)
{
  waveletFilter=wave.filter(wf)
  #Center ofEnergy of the scaling filter
  cesf=computeCenterEnergy(waveletFilter$length,waveletFilter$lpf);
  #Center ofEnergy of the wavelet filter
  cewf=computeCenterEnergy(waveletFilter$length,waveletFilter$hpf);
  return (c(cesf,cewf))
  
}

#auxiliar function: gets the center of energy of coefs. length(coefs)=l
computeCenterEnergy=function(l,coefs){
  auxVector=0:(l-1)
  coefs2=(coefs)^2
  a=coefs2%*%auxVector
  b=sum(coefs2)
  return (a/b)
}