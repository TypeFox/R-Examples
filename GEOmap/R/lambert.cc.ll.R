`lambert.cc.ll` <-
function(x,y, PROJ.DATA)
  {
###  lambert conformal conic Snyder(USGS) p. 104

    phi0 = PROJ.DATA$LAT0
    lam0 = PROJ.DATA$LON0
    phi1 = PROJ.DATA$LAT1
    phi2 = PROJ.DATA$LAT2
	FE = PROJ.DATA$FE
	FN = PROJ.DATA$FN
 
#  Constants:
    DR = pi/180
phi1=phi1*DR
phi2=phi2*DR
phi0=phi0*DR
lam0=lam0*DR
R = MAPconstants()$A.MAPK

x=x-FE
y=y-FN

n=log(cos(phi1)/cos(phi2))/log(tan(pi/4+phi2/2)/tan(pi/4+phi1/2))  #15-3
## print(paste(sep=' ', "n=", n))
F=cos(phi1)*(tan(pi/4+phi1/2))^n/n                                 #15-2
## print(paste(sep=' ', "F=", F))

rho0=R*F/(tan(pi/4+phi0/2))^n                                      #15-1a
##print(paste(sep=' ', "rho0=", rho0))

#  Calc rho,theta, phi, and lambda

rho=sign(n)*(x^2+(rho0-y)^2)^(1/2)                                 #14-10
## print(paste(sep=' ', "rho=", rho))


theta=atan(x/(rho0-y))                                             #14-11
## print(paste(sep=' ', "theta=", theta))

lam=theta/n+lam0                                                   #14-9
phi=2*atan((R*F/rho)^(1/n))-pi/2                                   #15-5

lon=(lam)*180/pi
lat=(phi)*180/pi
    
    return(list(lat=lat, lon=lon))
  }

