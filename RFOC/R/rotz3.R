`rotz3` <-
function( deg )
  {
    rad1 = deg * 0.0174532925199;
    r = diag(3)
    r[1, 1] = cos(rad1)
    r[1, 2] = sin(rad1)
    r[2, 2] = r[1, 1]
    r[2, 1] = -r[1, 2]
    return(r)

  }

