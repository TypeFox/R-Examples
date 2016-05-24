# Visualize P-splines
demo.PSplines <- function(y=NULL, x=NULL, ...)
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
# Simulate data
  if (is.null(y))
  {  
     n <- 100
     x <- seq(0, 1, length = n)
   set.seed(123)
     y <- 1.2 + sin(5  * x) + rnorm(n) * 0.2
    xg <- seq(0, 1, length = 500)
  xmin <- min(x)
  xmax <- max(x)
  }
  else 
  {
     y <- y
     x <- if (is.null(x)) stop ("the x-variable must be set here") else x
  xmin <- min(x)
  xmax <- max(x)
    xg <- seq(xmin, xmax, length = 500)  
  }
   ord <- 2
 order <- 2
   lla <- 1
  nseg <- 10
  bdeg <- 3
  pord <- 2
if (interactive()) 
    {
    ps.smooth = function(panel)
     {
      nseg <- floor(panel$nseg)
      bdeg <- panel$bdeg
       lla <- panel$lla
    lambda <- 10 ^ lla
      pord <- panel$ord
         B <- bbase(x,  xl = xmin, xr = xmax, nseg = nseg, deg = bdeg)
        nb <- ncol(B)
         D <- diff(diag(nb), diff = pord)
         P <- lambda * t(D) %*% D
         a <- solve(t(B) %*% B + P, t(B) %*% y)
         a <- as.vector(a)
        nb <- ncol(B)
      cols <- hcl(h = seq(60, 240, length = nb), c =90, l = 70)
        Bg <- bbase(xg, xl = xmin, xr = xmax, nseg = nseg, deg = bdeg)
         A <- diag(a)
         z <- Bg %*% a
         plot(x, y, ...) 
         matlines(xg, Bg %*% A, type = 'l', lty = 1, lwd = 2, col= cols, xlab = '', ylab = '', ylim = c(0, 1))
         lines(xg, z, col = 'red', lwd = 3)        
     knots <- seq(0, 1, length = nseg + 1)
         points(knots, 0 * knots, pch = 15, cex = 0.8)
        tl <- paste('P-splines, n = ', nb, ', order = ', pord, ', degree = ', bdeg, ', log10(lambda) = ', lla, sep = '')
        title(tl)
  panel
}

         ps.panel = rp.control('PSP', size = c(400, 200),  lla = 1, ord = 2, nseg = 10, bdeg = 3)
         rp.slider(ps.panel, variable = lla, from = -4, to = 6,  action = ps.smooth, resolution = 0.2, showvalue = T, title = 'Set log10(lambda)')
         rp.slider(ps.panel, variable = nseg, from = 3, to = 20, action = ps.smooth, resolution = 1, showvalue = T, title = 'Size of basis')
         rp.doublebutton(ps.panel, variable = ord,  action = ps.smooth, initval = 2,  step = 1, range = c(1, 4), showvalue = T, "Penalty order")
         rp.doublebutton(ps.panel, variable = bdeg, action = ps.smooth, initval = 3,  step = 1, range = c(0, 4), showvalue = T, "B-spline degree")
    }
}
