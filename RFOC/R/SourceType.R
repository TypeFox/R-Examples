SourceType <- function(v)
  {
### from Tape and Tape
###  given eigenvalues, get source type

###  step 2 in:
    ##   Miscellaneous notes for
    ##   Tape and Tape (2012a): “A geometric setting for moment tensors”
    ##    Carl Tape (carltape@gi.alaska.edu) October 19, 2012
                                        # http://www.giseis.alaska.edu/input/carl/research/beachball.html

    ##  equation 4 from
    ## Walter Tape and Carl Tape
    ## A geometric comparison of source-type plots for moment tensors
    ## Geophys. J. Int. (2012) 190, 499–510 doi: 10.1111/j.1365-246X.2012.05490.x

    v = sort(v, decreasing=TRUE)
    
    sq3 = sqrt(3)
    
    rho = sqrt(sum(v^2))


    vdif = v[1]-v[3]
    
    if(vdif==0)
      {
        lambda = pi
      }
    else
      {
        lambda = atan( (-v[1] +2*v[2] - v[3] )/(sq3*(vdif)))
      }

    cosj =  sum(v)/(sq3*rho)

    
    if(cosj>1 | cosj<(-1) )
      {
        
        if(cosj>1) { beta = 0 }
        if(cosj<(-1)) { beta = pi } 
      }
    else
      {
        beta = acos(cosj )
      }

    return(list(lam=lambda*180/pi, phi=beta*180/pi))
  }
