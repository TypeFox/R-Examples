cirrange = function( ang, radians=FALSE) {

                                        #  determine the additive constant.

  if(radians)
    cnst = pi * 2 
  else
    cnst = 360

                                        # deal with the lower limit.

  ang = ang %% cnst

                                        # deal with negative values, if any
  
  neg = which(ang < 0)
  ang[neg] = ang[neg] + cnst
  
  return(ang)
}
