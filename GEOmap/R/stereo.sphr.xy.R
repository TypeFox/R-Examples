`stereo.sphr.xy` <-
function(phi,lam, PROJ.DATA)
  {
	###  lambert conformal conic Snyder(USGS) p. 157

    lam = RPMG::fmod(lam, 360)	

    
    phi1 = PROJ.DATA$LAT0
    lam0 = PROJ.DATA$LON0
###    phi1 = PROJ.DATA$LAT1
    
    FE = PROJ.DATA$FE
    FN = PROJ.DATA$FN

    ##  print(paste(sep=' ', phi,lam, phi0, lam0,  phi1,  phi2, FE, FN))

                                        #  
    phi =phi*pi/180
    lam =lam*pi/180
    phi1=phi1*pi/180 
    lam0=lam0*pi/180

    R = MAPconstants()$A.MAPK

    k0 = 1
    k = 2*k0/(1+sin(phi1)*sin(phi) + cos(phi1)*cos(phi)*cos(lam-lam0))  ###  21-4

    x = R*k*cos(phi)*sin(lam-lam0)                                             #21-2
    y = R*k*(cos(phi1)*sin(phi)-sin(phi1)*cos(phi)*cos(lam-lam0))        #  21-2
    ##  print(paste(sep=' ',rho, theta, x, y))

    return(list(x=x, y=y))
  }

