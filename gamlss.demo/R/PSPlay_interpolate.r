# Visualize Whittaker smoother  with interpolations and extrapolation
demo.interpolateSmo <- function(y=NULL, w=NULL, ...)
{
  whitsm = function(y, lambda = 10, d = 2, w = 0 * y + 1){
    m <- length(y)
    E <- diag(m)
    D <- diff(E, diff = d)
    z <- solve(diag(w) + lambda * t(D) %*% D, w * y)
  return(z)
}
# Simulate data
  if (is.null(y))
  {  
     n <- 100
     x <- seq(0, 1, length = n)
    set.seed(123)
     y <- 1.2 + sin(6  * x) + rnorm(n) * 0.1
     w <- 0 * y + 1
     w[c(1:10, 45:65, 85:100)] <- 0
  }
  else 
  {
     y <- y
     w <- if (is.null(w)) stop ("the weights must be set here") else w 
     x <- seq(1, length(y), step=1)
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
        z <- whitsm(y, lambda, order, w= w)
       tl <- paste('Whittaker smoothing; order = ', order, ', log10(lambda) = ', lla, sep = '')
      sel <- w == 1         
     plot(x[sel], y[sel], pch = 15, cex = 0.8, xlab = '', ylab = '', main = tl, ...)
     lines(x, z, col = 'red', lty = 1, lwd = 3)
      panel
     }
         ws.panel = rp.control('Whittaker', lla = 1)
         rp.slider(ws.panel, variable = lla, from = -2, to = 10, action = ws.smooth, resolution = 0.1, showvalue = T, title = 'Set log10(lambda)')
         rp.doublebutton(ws.panel, variable = order, action = ws.smooth, initval = 2, step = 1, range = c(1, 4), showvalue = T, "Order")
    }
}
