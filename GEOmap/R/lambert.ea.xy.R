`lambert.ea.xy` <-
function(phi,lam, PROJ.DATA)
  {
	###  lambert azimuthal equal area Snyder(USGS) p. 185

lam = RPMG::fmod(lam, 360)	

 
    phi0 = PROJ.DATA$LAT0
    lam0 = PROJ.DATA$LON0
    phi1 = PROJ.DATA$LAT1
    phi2 = PROJ.DATA$LAT2
	FE = PROJ.DATA$FE
	FN = PROJ.DATA$FN

##  print(paste(sep=' ', phi,lam, phi0, lam0,  phi1,  phi2, FE, FN))

#  Constants:
phi =phi*pi/180
lam =lam*pi/180
phi1=phi1*pi/180 
phi2=phi2*pi/180 
phi0=phi0*pi/180
lam0=lam0*pi/180
R = MAPconstants()$A.MAPK
##  Az = lam
##  c = 

## print(paste(sep=' ',R, n, F, rho0))


## rho = 2*R*sin(c/2)                       # 24-1
##  theta = pi-Az                            # 20-2
##  hp = cos(c/2)                            # 24-1a
##  kp = sec(c/2)                            # 24-1b


kp = sqrt(2/(1 + sin(phi1)*sin(phi)+cos(phi1)*cos(phi)*cos(lam-lam0)))  ## 22-2

x = R*kp*cos(phi)*sin(lam-lam0) +FE         #22-4



y = R*kp*(cos(phi1)*sin(phi)-sin(phi1)*cos(phi)*cos(lam-lam0))  + FN ## 22-5
##  print(paste(sep=' ',rho, theta, x, y))

    return(list(x=x, y=y))
  }

