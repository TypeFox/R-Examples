`gclc` <-
function(phiorg, lamorg, phi, lam)
{
 
 A=6378206.4
 E2=0.00676866
 E=0.0822719
 E1=0.993231340
 TwoE=0.164543800
 R=6378.2064
 DEG2RAD= 0.017453293
 RAD2DEG= 57.295778667
 EARTHRAD= 6378.163
 ECCEN= 0.0033670033
 PO180= 0.017453293

 
  PHI0= DEG2RAD * phiorg
  LAM0= DEG2RAD * lamorg
 
  crlat= atan(ECCEN * sin(2.0 * PHI0) / GCLCFR(PHI0))
  zgcl= PHI0 - crlat
  a= GCLCFR(zgcl)
  rho= EARTHRAD * a
  b= ((ECCEN * sin(2.0 * zgcl) / a))^2 + 1.0
  ca= 2.0 * ECCEN * cos(2.0 * zgcl) * a + ((ECCEN * sin(2.0 * zgcl)))^2
  cdist= ca / (a*a * b) + 1.0
  RHO0= rho
  C= crlat
  B= cdist
  
  xlat= DEG2RAD * phi
  ylon= DEG2RAD * lam
  y= (-1.0)*RHO0 * (LAM0 - ylon) * cos(xlat - C)
  x= RHO0 * (xlat - PHI0) / B
  
list(x=y,y=x)


}

