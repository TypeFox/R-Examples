#////////////////////////////
#// Seed dispersal function (Colbach and al. 2000b)
#////////////////////////////
# Example:
# distance = seq(1,1.5, by=0.05)
# plot(x=distance, y = fseed(distance))

fseed=function  ( point) 
{
 distance= point[1];
 grainesB =  1.38930;
 grainesC =  2.08686;

   fd = (grainesB * grainesC) * distance **(grainesC-2) *
         exp(-grainesB * distance ** (grainesC)) / (2.0*pi);
 
  return(fd);
}

#////////////////////////////
#// f1: Pollen dispersal function (Klein)
#////////////////////////////
# Example:
#a=matrix(distance, ncol=1)
#b= apply(a,1,fpollen)
# plot(x=distance, y =b)
# lines(x=distance, y = fseed(distance))

fpollen = function ( point)
{
 r = point[1];

  if ( r<=1.5 ) {
    rt = 0.340-0.405*r+0.128*r**(2.0);
  }  else {
    if (r <= 50) {
      rt = 0.03985/(1+r**(3.12)/3.80);
    }
    else {
      rr = 50
    gamma = -2.29
#    //heaviest value: -2.14, lighest: - 2.56
    K = (0.03985 / (1.0+rr**(3.12)/3.80)) / (1.0+ rr)**gamma
    rt = K * (1.0 + r)**gamma
}      
  }

  return (rt);
} 
#////////////////////////////
#// Function constant
#////////////////////////////
fcte = function ( distance)
{
return(1)
}
