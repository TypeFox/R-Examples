JuliaImage <-
function(imageN,centre,L,C) {
#  Generates Julia Set on the given rectangle 
#  on the complex plane.
#
#       z <- z^2 + C 
#
# Returns image matrix
# 
#   imageN : pixel size NxN
#   centre : centre coordinate of the image
#   L      : side length 
#   C      : constant to be used to generate Julia set
#
   delta <- L/(imageN)
   image <- array(0,c(imageN,imageN));
   Re <- Re(centre) - L/2.0 ;# Re is on the x-axis
   for(i in 1:imageN) {
      Im <- Im(centre) + L/2.0 ;# Im is on the y-axis
     for(j in 1:imageN) {
        z <- (Re + Im*1i)
        image[j,i] <- JuliaIterate(z,C)
        Im <- Im - delta 
      }
        Re <- Re + delta 
   }
   image <- image/max(image);# scale it. 
 return(image);
}

