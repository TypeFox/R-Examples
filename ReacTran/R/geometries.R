## =============================================================================
##
## Surface area and volume of geometries
##
## =============================================================================


# surface and volume of a sphere
  g.sphere <- function (x) return(list(surf=4*pi*x*x,vol=4/3*pi*x^3))
 
# surface and volume of a spheroid
  g.spheroid <- function (x, b=1) {     # b = long/short radius
    bb <- x*b
    te <- acos(b)
    if (b < 1)      #oblate
      surf <- 2*pi *(x^2 + bb^2 / sin(te)*log((1+sin(te))/cos(te)))
    else if (b > 1)
      surf <- 2*pi*(x^2 + x*bb*te/sin(te))
    else
      surf <- 4 * pi*x^2

   vol <- 4/3*pi * x^2 * bb
   return(list(surf=surf,vol=vol))
  }

# surface and volume of a cylinder 
  g.cylinder <- function (x,L=1) return(list(surf=2*pi*x*L,vol=pi*x^2*L))
