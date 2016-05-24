`stereo.sphr.ll` <-
function(x,y, PROJ.DATA)
  {
###  lambert conformal conic Snyder(USGS) p. 104
    k0 = 1
    
    lam0 = PROJ.DATA$LON0
###    phi1 = PROJ.DATA$LAT1
    phi1 = PROJ.DATA$LAT0
    
                                        #  Constants:
    phi1=phi1*pi/180 
    lam0=lam0*pi/180
    R = MAPconstants()$A.MAPK



    rho =sqrt(x^2+y^2)                                       #20-18

                                        #  Calc rho,theta, phi, and lambda

    cee = 2*atan2( rho, 2*R*k0 )                               #21-15
    ## print(paste(sep=' ', "rho=", rho))


    phi=asin(cos(cee)*sin(phi1) +y*sin(cee)*cos(phi1)/rho)                                #20-14

    if(PROJ.DATA$LAT1==90)
      {
        lam=lam0+atan2(x,(-y))  #20-15
      }
      if(PROJ.DATA$LAT1==(-90) )  
      {
        lam=lam0+atan2(x,(y))  #20-16
        
      }
    if(PROJ.DATA$LAT1!=(-90) & PROJ.DATA$LAT1!=(90) )
      {
        lam=lam0+atan2(x*sin(cee),   rho*cos(phi1)*cos(cee)-y*sin(phi1)*sin(cee) )   ### 20-15
        
      }
    
    
    
    
    lon=(lam)*180/pi
    lat=(phi)*180/pi
    
    return(list(lat=lat, lon=lon))
  }

