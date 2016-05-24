## =============================================================================
## Creates the Runge-Kutta figure from the initial value problem chapter
## =============================================================================

# Take an explicit Runge-Kutta step and return the coefficients 
rkfixed <- function(y, t, h, f, A, b1, b2, c, k) {            

  k[1] <- f(t, y)
  n  <- length (k)
  yn <- NULL 
  xn <- NULL  
  for (i in 1:n) 
      k[i] <- f (t + c[i]*h , y + h * sum(A[i,1:(i-1)]*k[1:(i-1)]))
  for (i in 1:n) {
      yn <- c(yn, y + h * sum(A[i,1:(i-1)]*k[1:(i-1)]))
      xn <- c(xn, t + c[i]*h)
  }
  ynext <- y + h * sum(b1 * k)
  error <- max( abs(h * sum((b2-b1) * k)))
  list (y = ynext, yerr = error, k = k, yn = yn, xn=xn)                                               
}

# Take an implicit Runge-Kutta step
rkimp <- function(y, t, h, f, A, b1, b2, c, k) {            
  n     <- length (k)
  kdiff <- numeric(n)
  
  rootfun <- function(kn) {
    for (i in 1:n) 
      kdiff[i]<-  -kn[i] + f (t + c[i]*h , y + h * sum(A[i,]*kn) )
    return (kdiff)
  }
  xn <- NULL
  yn <- NULL
  k     <- multiroot(f=rootfun, start=runif(n))$root
  for (i in 1:n)  {
    yn <- c(yn, y + h * sum(A[i,]*k))
    xn <- c(xn, t + c[i]*h)
  }
  ynext <- y + h * sum(b1 * k)
  list (y = ynext, k=k, yn = yn, xn = xn)                                               
}

# Function to decide if explicit or implicit method
FullY <- function(y, t, h, f, A, b1, b2, c, k) {
    if (sum(A [upper.tri(A)]) !=  0)  # is implicit method
      RK <- rkimp(y, t, h, f, A, b1, b2, c, k)
    else 
      RK <- rkfixed(y, t, h, f, A, b1, b2, c, k)
    list(xn = RK$xn, yn = RK$yn, k = RK$k, xx = t+h, yy = RK$y)
  }

###################################################
### Function to plot rungekuttas
###################################################

rkMethodPlot <-  function (rk, ...) {

  A    <- rk$A
  b1   <- rk$b1
  b2   <- rk$b2
  c    <- rk$c

  dots   <- list(...)
  nmdots <- names(dots)

  if(is.null(dots$main)) dots$main <- rk$ID
  i1 <- ifelse (c[1] == 0, 0, 1 ) 

  n <- nrow(A)
  k <- vector(length=n)

  yini <- 0.01
  mu   <- 0.1
  # toy derivative function
  f <- function (t, y) mu*y

  par(mar=c(3,4,3,2))
  x1  <- 1
  y1  <- yini*exp(mu*x1)
  y2  <- yini*exp(mu*40)

  h   <- 10
  YY  <- FullY(yini, 0, h, f, A, b1, b2, c, k)

  do.call("curve", c(alist(yini * exp(mu*x), -1, 11, lwd = 1,xlab = "",
        ylab = "", axes = FALSE, frame.plot = TRUE,
        ylim = c(0.005, 0.03)), dots))


  axis(side = 2, at = yini*exp(mu*c(0,10)), 
       labels =c(expression(y[0]), expression(y[1])))
  axis(side = 1, at = c(0,10), 
       labels =c(expression(x[0]), expression(x[1])))

  points(YY$xn, YY$yn)
  points(YY$xx, YY$yy)
  points(0, yini, pch = 16, cex = 1.5, col = "grey")
  points(YY$xx,YY$yy, pch = 16, cex = 1.5, col = "black")
  if (i1) 
    text( YY$xn, YY$yn-0.001, 1:n)
  else
    text( YY$xn, YY$yn-0.001, 0:(n-1))

  dx <- 0.5
  x  <- YY$xn
  yy <- YY$yn
  y1 <- yy-dx*YY$k
  y2 <- yy+dx*YY$k
  segments(x-dx, y1, x+dx, y2, lwd=2)
  return(YY)
}
