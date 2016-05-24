inverseTAPE<-function(GAMMA, BETA)
  {
####  angles in degrees
###  GAMMA = long
###  BETA  = lat 
    ##   Inverse Tape plot from Moment tensor
    sq3 = sqrt(3)
    G = sq3*tan(GAMMA*pi/180)
    B = sq3*cos(BETA*pi/180)

    G2 = G^2
    B2 = B^2

    G23 = 3+G2

    k2 = (-12*B - 4*B*G2)^2 - 4*(18+6*G2)*(-9+5*B2 + 6*G - 2*B2*G - G2 +B2*G2)

    if(k2>0)
      {
        bigrad =
          sqrt( k2  )
      }
    else
      {
        bigrad = 0
      }

    
    LAM1a = (12*B+4*B*G2 - bigrad )/(12*G23)

    LAM2a = (B - B*G + ((2*B*G)/(G23) + (2*B*G^3)/(3*(G23)) - G* bigrad/(6*(G23))))/(3-G)

    LAM3a = (-2*B + 3*B/(G23) +(B*G)/G23 +(B*G2)/G23+ (B*G^3)/(3*G23) - bigrad/(4*G23) - (1/(12*G23)) * G*bigrad )/(G-3)

    Va = c(LAM1a, LAM2a, LAM3a)
    Va = Va/max(abs(Va))

    

    LAM1b = (12*B+4*B*G2 + bigrad )/(12*G23)

    LAM2b = (B - B*G + ((2*B*G)/(G23) + (2*B*G^3)/(3*(G23)) + G* bigrad/(6*(G23))))/(3-G)

    LAM3b = (-2*B + 3*B/(G23) +(B*G)/G23 +(B*G2)/G23+ (B*G^3)/(3*G23) + bigrad/(4*G23) + (1/(12*G23)) * G*bigrad )/(G-3)

    
    Vb = c(LAM1b, LAM2b, LAM3b)
    Vb = Vb/max(abs(Vb))

    return(list(Va=Va, Vb=Vb) )

  }


