#////////////////////////////
#// fonction de dispersion des graines (Colbach et al. 2000b)
#////////////////////////////
# Exemple:
# distance = seq(1,1.5, by=0.05)
# plot(x=distance, y = fseed(distance))

fseed=function  ( distance) 
{
 grainesB =  1.38930;
 grainesC =  2.08686;

   fd = (grainesB * grainesC) * distance **(grainesC-2) *
         exp(-grainesB * distance ** (grainesC)) / (2.0*pi);
 
  return(fd);
}

#////////////////////////////
#// f1: la fonction de dispersion du pollen de Klein
#////////////////////////////
# Exemple:
#a=matrix(distance, ncol=1)
#b= apply(a,1,fpollen)
# plot(x=distance, y =b)
# lines(x=distance, y = fseed(distance))

fpollen = function ( distance)
{
 r = distance

  if ( r<=1.5 ) {
    rt = 0.340-0.405*r+0.128*r**(2.0);
  }  else {
    if (r <= 50) {
      rt = 0.03985/(1+r**(3.12)/3.80);
    }
    else {
      rr = 50
    gamma = -2.29
#    //valeur la plus lourde -2.14, + légère - 2.56
    K = (0.03985 / (1.0+rr**(3.12)/3.80)) / (1.0+ rr)**gamma
    rt = K * (1.0 + r)**gamma
}      
  }

  return (rt);
} 
