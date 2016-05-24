demo.BSplines <- function()
{
 tpower <- function(x, t, p)
 # Truncated p-th power function
    (x - t) ^ p * (x > t)
#--------------
  bbase <- function(x, xl = min(x), xr = max(x), nseg = 10, deg = 3)
  {
    # Construct B-spline basis
       dx <- (xr - xl) / nseg
    knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
        P <- outer(x, knots, tpower, deg)
        n <- dim(P)[2]
        D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
        B <- (-1) ^ (deg + 1) * P %*% t(D)
        B
  }
#--------------
    nseg <- 5
    bdeg <- 3
     x <- seq(0, 1, length = 500)
 order <- 2
   lla <- 1
if (interactive()) 
    {
    #---------------
 bs.random = function(panel) 
 {
     nseg <- floor(panel$nseg)
     bdeg <- panel$bdeg
        B <- bbase(x, nseg = nseg, deg = bdeg)
        a <- runif(ncol(B))
       nb <- ncol(B)
     cols <- hcl(h = seq(60, 240, length = nb), c =90, l = 70)
        A <- diag(a)
       z <- B %*% a
  matplot(x, B %*% A, type = 'l', lty = 1, lwd = 2, col= cols, xlab = '', ylab = '', ylim = c(0, 1))
  lines(x, z, col = 'red', lwd = 3)        
  knots <- seq(0, 1, length = nseg + 1)
  points(knots, 0 * knots, pch = 15, cex = 0.8)
  title(paste('B-spline basis, n = ', nb, ', degree = ', bdeg, sep = ''))
  return(panel)
 }
#---------------
    bs.draw = function(panel)
    {
     nseg <- floor(panel$nseg)
     bdeg <- panel$bdeg
        B <- bbase(x, nseg = nseg, deg = bdeg)
        a <- rep(1, ncol(B))
       nb <- ncol(B)
     cols <- hcl(h = seq(60, 240, length = nb), c =90, l = 70)
        A <- diag(a)
        z <- B %*% a
  matplot(x, B %*% A, type = 'l', lty = 1, lwd = 2, col= cols, xlab = '', ylab = '', ylim = c(0, 1))
  lines(x, z, col = 'red', lwd = 3)        
    knots <- seq(0, 1, length = nseg + 1)
  points(knots, 0 * knots, pch = 15, cex = 0.8)
  title(paste('B-spline basis, n = ', nb, ', degree = ', bdeg, sep = ''))
      panel
    }
         bs.panel <- rp.control('B-spline visualizer', nseg = 5, size = c(300, 100))
         rp.slider(bs.panel, nseg, 1, 20, action=bs.draw, "B-splines")
         rp.doublebutton(bs.panel, variable = bdeg, action = bs.draw, initval = 3, step = 1, range = c(0, 4), showvalue = T, "Degree")
         rp.button(bs.panel, action = bs.random, 'Random')     
    }
   
}
