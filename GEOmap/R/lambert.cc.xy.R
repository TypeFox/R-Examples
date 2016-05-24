`lambert.cc.xy` <-
function(phi,lam, PROJ.DATA)
  {
	###  lambert conformal conic Snyder(USGS) p. 104

lam = RPMG::fmod(lam, 360)	

 
    phi0 = PROJ.DATA$LAT0
    lam0 = PROJ.DATA$LON0
    phi1 = PROJ.DATA$LAT1
    phi2 = PROJ.DATA$LAT2
	FE = PROJ.DATA$FE
	FN = PROJ.DATA$FN

##  print(paste(sep=' ', phi,lam, phi0, lam0,  phi1,  phi2, FE, FN))

#  Constants:

DR = pi/180
phi =phi*DR
lam =lam*DR
phi1=phi1*DR
phi2=phi2*DR
phi0=phi0*DR
lam0=lam0*DR
R = MAPconstants()$A.MAPK

n=log(cos(phi1)/cos(phi2))/log(tan(pi/4+phi2/2)/tan(pi/4+phi1/2))  #15-3
F=cos(phi1)*(tan(pi/4+phi1/2))^n/n                                 #15-2
rho0=R*F/(tan(pi/4+phi0/2))^n                                      #15-1a

## print(paste(sep=' ',R, n, F, rho0))


rho = R*F/((tan(pi/4+phi/2))^n)                                    #15-1
theta = n*(lam-lam0)                                               #14-4
x = rho*sin(theta)+FE                                              #14-1
y = rho0-rho*cos(theta)+FN                                         #14-2
##  print(paste(sep=' ',rho, theta, x, y))

    return(list(x=x, y=y))
  }

