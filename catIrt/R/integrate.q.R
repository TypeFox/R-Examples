# Quick integration:
integrate.q <-
function(x, y){

  n <- length(x)
  
# Shift to find lower and upper vectors: 1-(n - 1) and 2-n
  lx <- x[-n]; ux <- x[-1]
  ly <- y[-n]; uy <- y[-1]
 
  sum( (ux - lx) * (ly + uy) / 2 )

} # END integrate.q