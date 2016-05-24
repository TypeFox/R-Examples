co_refract_forward=function( a, p=1010, t=283) {
  d2r = pi/180.

  w = (a<15)
  r = 0.0166667/tan((a + 7.31/(a+4.4))*d2r)
  r[w] = 3.569*(0.1594 + .0196*a[w] +
     .00002*a[w]^2)/(1.+.505*a[w]+.0845*a[w]^2)
  tpcor = p/1010. * 283/t
  r = tpcor * r
  return(r)
}