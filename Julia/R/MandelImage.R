MandelImage <-
function(imageN,centre,L) {
#  Generates Mandelbrot Set on the given rectangle
#  on the complex plane.
# 
#       z <- z^2 + z_0 
#
# Returns image matrix
# 
#   imageN : pixel size NxN
#   centre : centre coordinate of the image
#   L      : side length 
#
   delta <- L/(imageN)
   image <- array(0,c(imageN,imageN));
     Re <- Re(centre) - L/2.0 ;# Re is on the x-axis
   for(i in 1:imageN) {
     Im <- Im(centre) + L/2.0 ;# Im is on the y-axis
    for(j in 1:imageN) {
        z_0 <- (Re + Im*1i)
        image[j,i] <- MandelIterate(z_0)
        Im <- Im - delta 
      }
       Re <- Re + delta 
    }
  image <- image/max(image)
  return(image)
}

