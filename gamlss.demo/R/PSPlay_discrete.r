# Visualize Whittaker smoother
demo.RandomWalk <- function(y=NULL,...)
{
  whitsm = function(y, lambda = 10, d = 2)
     {
     m <- length(y)
     E <- diag(m)
     D <- diff(E, diff = d)
     z <- solve(E + lambda * t(D) %*% D, y)
     return(z)
     }
# Simulate data
  if (is.null(y))
  {  
     n <- 100
     x <- seq(0, 1, length = n)
     set.seed(123)
     y <- 1.2 + sin(5  * x) + rnorm(n) * 0.2
  }
  else 
  {
     y <- y
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
        z <- whitsm(y, lambda, order)# gamlss(y~rw(1))
       tl <- paste('Whittaker smoothing; order = ', order, ', log10(lambda) = ', lla, sep = '')
       plot(x, y, type = 'l', xlab = '', ylab = '', lwd = 1.5, main = tl, ...)
     lines(x, z, col = 'red', lty = 1, lwd = 3)
      panel
   }
         ws.panel = rp.control('Whittaker', lla = 1)
         rp.slider(ws.panel, variable = lla, from = -2, to = 10, action = ws.smooth, resolution = 0.1, showvalue = T, title = 'Set log10(lambda)')
         rp.doublebutton(ws.panel, variable = order, action = ws.smooth, initval = 2, step = 1, range = c(1, 4), showvalue = T, "Order")
  
   }
}
