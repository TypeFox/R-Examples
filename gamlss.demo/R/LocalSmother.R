# Visualize local mean smoother
demo.Locmean <- function(y=NULL, x=NULL, ...)
{
#--------------
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
     x <- if (is.null(x)) stop ("the x-variable must be set here") else x
  }
 span  <- 0.5
if (interactive()) 
    {
    ps.smooth = function(panel)
     {
     span <- panel$span
         plot(x, y, ...) 
         m1 <- Locmean(y,x, span=span)
         lines(x, fitted(m1), col = 'red', lwd = 3)        
        tl <- paste('Local mean smoothing , span = ', span, sep = '')
        title(tl)
  panel
}
         ps.panel = rp.control('Local Mean Smoother', size = c(400, 200),  span = 0.5)
         rp.slider(ps.panel, variable = span, from = 0.01, to = 2,  action = ps.smooth, resolution = 0.01, showvalue = T, title = 'span')
    }
}
#---------------------------------------------------------------
demo.Locpoly <- function(y=NULL, x=NULL, ...)
{
#--------------
# Simulate data
  if (is.null(y))
  {  
    n <- 100
    x <- seq(0, 1, length = n)
   set.seed(123)
    y <- 1.2 + sin(5  * x) + rnorm(n) * 0.2
   #xg <- seq(0, 1, length = 500)
  }
  else 
  {
     y <- y
     x <- if (is.null(x)) stop ("the x-variable must be set here") else x
 # xmin <- min(x)
 # xmax <- max(x)
 #   xg <- seq(xmin, xmax, length = 500)  
  }
 span  <- 0.5
 order <- 2
if (interactive()) 
    {
    ps.smooth = function(panel)
     {
     span <- panel$span
     order <-  panel$order
         plot(x, y, ...) 
         m1 <- Locpoly(y,x, span=span, order=order)
         lines(x, fitted(m1), col = 'red', lwd = 3)        
        tl <- paste('Local Poly smoothing , span = ', span, sep = '')
        title(tl)
  panel
}

         ps.panel = rp.control('Local Polynomial Smoother', size = c(400, 200), span=0.5)
         rp.slider(ps.panel, variable = span, from = 0.1, to = 2,  action = ps.smooth, resolution = 0.1, showvalue = T, title = 'span')
         #rp.doublebutton(ps.panel, variable = span,  action = ps.smooth, initval = 0.5,  step = 0.1, range = c(0, 2), showvalue = T, "Polynomial order")
         rp.doublebutton(ps.panel, variable = order,  action = ps.smooth, initval = 2,  step = 1, range = c(1, 4), showvalue = T, "Polynomial order")
    }
}
#---------------------------------------------------------------
demo.WLocmean <- function(y=NULL, x=NULL, ...)
{
#--------------
# Simulate data
  if (is.null(y))
  {  
    n <- 100
    x <- seq(0, 1, length = n)
   set.seed(123)
    y <- 1.2 + sin(5  * x) + rnorm(n) * 0.2
 #  xg <- seq(0, 1, length = 500)
  }
  else 
  {
     y <- y
     x <- if (is.null(x)) stop ("the x-variable must be set here") else x
  #xmin <- min(x)
  #xmax <- max(x)
  #  xg <- seq(xmin, xmax, length = 500)  
  }
 lambda  <- 0.5
if (interactive()) 
    {
    ps.smooth = function(panel)
     {
         lambda <- exp(panel$lambda)
         plot(x, y,  ...) 
             m1 <- WLocmean(y, x, lambda=lambda)
            lines(x, fitted(m1), col = 'red', lwd = 3)        
            tl <- paste('Kernel smoothing , lambda = ', signif(lambda, digits = 2), sep = '')
        title(tl)
  panel
}
         ps.panel = rp.control('Locally Weighed Mean Smoother', size = c(400, 200), lambda=.5)
         rp.slider(ps.panel, variable = lambda, from = -6, to = 6,  action = ps.smooth, resolution = 0.2, showvalue = T, title = 'Set log(lambda)')
        # rp.doublebutton(ps.panel, variable = lambda,  action = ps.smooth, initval = 0.5,  step = 0.1, range = c(0, 2), showvalue = T, "lamba")
    }
}
#----------------------------------------------------------------------
#---------------------------------------------------------------
demo.WLocpoly <- function(y=NULL, x=NULL, ...)
{
#--------------
# Simulate data
  if (is.null(y))
  {  
    n <- 100
    x <- seq(0, 1, length = n)
   set.seed(123)
    y <- 1.2 + sin(5  * x) + rnorm(n) * 0.2
   #xg <- seq(0, 1, length = 500)
  }
  else 
  {
     y <- y
     x <- if (is.null(x)) stop ("the x-variable must be set here") else x
  #xmin <- min(x)
  #xmax <- max(x)
    #xg <- seq(xmin, xmax, length = 500)  
  }#
 lambda  <- 0.5
   order <- 2
if (interactive()) 
    {
    ps.smooth = function(panel)
     {
     lambda <- exp(panel$lambda)
      order <-  panel$order
         plot(x, y, ...) 
         m1 <- WLocpoly(y, x, lambda=lambda, order=order)
         lines(x, fitted(m1), col = 'red', lwd = 3)        
        tl <- paste('Local Poly smoothing , lambda = ', signif(lambda, digits = 2), sep = '')
        title(tl)
  panel
}

         ps.panel = rp.control('Local Weighted Polynomial Smoother', size = c(400, 200), , lambda=.5)
         rp.slider(ps.panel, variable = lambda, from = -6, to = 6,  action = ps.smooth, resolution = 0.2, showvalue = T, title = 'Set log(lambda)')
         rp.doublebutton(ps.panel, variable = order,  action = ps.smooth, initval = 2,  step = 1, range = c(1, 4), showvalue = T, "Polynomial order")
    }
}
#---------------------------------------------------------
