eqpole = function(l,b,southpole=F) {
  radeg = 180/pi
  if(southpole){
    l1 = -l/radeg
    b1 = -b/radeg
  } else {
    l1 = l/radeg
    b1 = b/radeg
  }
  
  sq = 2*(1 - sin(b1))

  sq[sq<0] = 0.0
  r = 18.0e0*3.53553391*sqrt(sq)
  y =r*cos(l1)
  x =r*sin(l1)

  return(list(x=x,y=y))
}
