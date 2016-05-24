`FixDip` <-
function(A)
  {
                                        # given an azimuth and an angle
    az  = A$az
    dip = A$dip
    DEG2RAD = pi/180;
    RAD2DEG = 180/pi;
    tdip = DEG2RAD*RPMG::fmod(dip, 360.);
    co = cos(tdip)
    si = sin(tdip)
    ang = RAD2DEG*atan2(si, co)


    quad = rep(1, length(dip))

    quad[co>=0 & si>=0] = 1
    quad[co<0 & si>=0] = 2
    quad[co<0 & si<0] = 3
    quad[co>=0 & si<0] = 4

    dip[quad==1] = ang[quad==1]
    dip[quad==2] = 180-ang[quad==2]
    dip[quad==3] = 180+ang[quad==3]
    dip[quad==4] = -ang[quad==4]



    az[quad==1] = az[quad==1]
    az[quad==2] = 180+az[quad==2]
    az[quad==3] = az[quad==3]
    az[quad==4] = 180+az[quad==4]
    

    A$az = RPMG::fmod(az, 360.);

    
    A$dip = dip
    
    return(A);
  }

