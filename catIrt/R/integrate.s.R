# Sympson's integration rules:
# S1 = (b - a)/6 * [f(a) + 4f((a + b)2) + f(b)]
# S2 = (b - a)/8 * [f(a) + 3f((2a + b)/3) + 3f((a + 2b)/3) + f(b)]

integrate.s <-
function(f, min = -4, max = 4, quad = 33, method = "S1", ...){

# Calculate base, eval, and length:
  x <- seq(min, max, length.out = quad)
  y <- f(x, ...)
  n <- quad
  
# Shift to find lower and upper vectors: 1-(n - 1) and 2-n
  lx <- x[-n]; ux <- x[-1]
  ly <- y[-n]; uy <- y[-1]
 
  if( method == "S1" )
    sum( (ux - lx) / 6 * ( ly + uy + 4 * f( (lx + ux)/2, ... ) ) )
  else
    sum( (ux - lx) / 8 * ( ly + uy + 3 * ( f( (2*lx + ux)/3, ... ) + f( (2*ux + lx)/3, ... ) ) ) )

} # END integrate.s