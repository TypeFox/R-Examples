`grotseis` <-
function(ang, flip=FALSE)
{
  if(missing(flip)) { flip = FALSE }
                                        # angle is in degrees
  a=ang*pi/180
  if(flip)
    {
      f =-1;
    }else
  {
    f = 1;
  }
                                        #  the minus 1 on the vertical is for flipping it.
  rot=matrix(c(f, 0, 0, 0, cos(a), sin(a), 0, -sin(a), cos(a)), ncol=3)

  return(rot)

}

