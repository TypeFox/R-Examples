# Visualize Whittaker smoother for histogram
#---------------------------------------------
# The Poisson smoother
demo.histSmo <- function(y=NULL,...)
{
  histsm <- function(y, lambda = 10, d = 2)
 {
  # Penalty stuff
  m <- length(y)
  E <- diag(m)
  D <- diff(E, diff = d)            
  P <- lambda * t(D) %*% D  
  # Initialie
  u <- log(y + 0.5)
 z0 <- 0
  # Iterate
  for (it in 1:20) 
    {
     z <- exp(u)
    dz <- max(abs(z - z0))
    z0 <- z
     Z <- diag(z)
    u  <- solve(Z + P, y - z + z * u)
    if (dz < 1e-5) break
    }
  return(z)
 }
# Simulate data
  if (is.null(y))
  {  
          n <- 100
          set.seed(123)
          v <- rnorm(n)
        hst <- hist(v, breaks = seq(-3, 3, by = 0.1), ylab="y")
  hst$xname <- "y"
          x <- hst$mids
          y <- hst$counts
  }
  else 
  {
      y <- y
     Ry <-range(y)
     Ey <-(Ry[2]-Ry[1])*.10
    hst <- hist(y, breaks = seq((Ry[1]-Ey), (Ry[2]+Ey), length = floor(length(y)/10)))
      x <- hst$mids
      y <- hst$counts
  }
 order <- 2
   lla <- 1
if (interactive()) 
    {
    ws.smooth = function(panel)
    {
      lla <- panel$lla
    order <- panel$order
   lambda <- 10 ^ lla
        z <- histsm(y, lambda, order)
       tl <- paste('Histogram smoothing; order = ', order, ', log10(lambda) = ', lla, sep = '')
       plot(hst, ...)
    lines(x, z, col = 'red', lty = 1, lwd = 3)
      panel
   }
         ws.panel = rp.control('Smooth Histogram', lla = 1)
         rp.slider(ws.panel, variable = lla, from = -2, to = 10, action = ws.smooth, resolution = 0.1, showvalue = T, title = 'Set log10(lambda)')
         rp.doublebutton(ws.panel, variable = order, action = ws.smooth, initval = 2, step = 1, range = c(1, 4), showvalue = T, "Order")
   }
}
