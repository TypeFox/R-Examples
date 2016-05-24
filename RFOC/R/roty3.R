`roty3` <-
function( deg )
  {
    rad1 = deg * 0.0174532925199;
    r = diag(3)
    r[1, 1] = cos(rad1)
    r[3, 1] = sin(rad1)
    r[3, 3] = r[1, 1]
    r[1, 3] = -r[3, 1]
    return(r)

  }

