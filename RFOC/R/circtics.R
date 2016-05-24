`circtics` <-
function(r=1, dr=.02, dang=10,... )
  {
    if(missing(r)) { r=1 }
    if(missing(dr)) { dr=0.02 }
    if(missing(dang)) { dang=10  }
 
    
    
# %    plot a set of tic marks and the center of the RPMG::circle
    #########

    
    points(0,0, pch=3)
    DR = r+dr

    ang =pi*seq(0, 360, by=dang)/180

    segments(r*cos(ang), r*sin(ang), DR*cos(ang), DR*sin(ang), ...)
    
}

