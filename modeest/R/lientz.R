# Author: P. Poncet

lientz <-
function(x,         # sample (the data)
         bw = NULL) # corresponds to parameter 'beta' in Lientz (1970)
{
#############################
# Empirical Lientz's function
#############################

  if (bw <= 0 | bw >= 1) stop("argument 'bw' must belong to (0, 1)")
  
  y <- sort(x)
  ny <- length(y)
  k <- ceiling(bw*ny) - 1
  
  if (k==0) {
    f <-
    function(z)
    {
      yy <- (y + c(y[-1],Inf))/2
      i <- sapply(z, FUN = function(zz) min(which(zz <= yy)))      
      return(abs(y[i] - z))    
    }
  
  } else if (k>0) {
    f <-
    function(z)
    {
      yy1 <- c(y[(k+1):ny],rep(Inf,k))   # k+1 lags
      yy2 <- c(y[(k+2):ny],rep(Inf,k+1))   # k+2 lags
      yy <- sort(c((y+yy1)/2,(y+yy2)/2))
      
      i <- sapply(z, FUN = function(zz) min(which(zz <= yy)))
      j <- i%%2
      yy <- y[j*k + (i+j)/2]
      return(ifelse(j==1, yy-z, z-yy))
    }
  }
  
  class(f) <- c("lientz", class(f))
  attr(f, "call") <- sys.call()
  attr(f, "x") <- x
  attr(f, "bw") <- bw
  attr(f, "source") <- NULL
  
  ## Output
  return(f)
  
}


#Lientz <- lientz


plot.lientz <-
function(x,            # an object of class 'lientz'
         zoom = FALSE, # if TRUE, one can zoom on the graph created
         ...)
{
  if (!inherits(x, "lientz")) {
    stop("argument 'x' must inherit from class 'lientz'")
  }

  arg <- list(...)
  ylim <- arg$ylim
  main <- arg$main
  xlab <- arg$xlab
  ylab <- arg$ylab


  xx <- attr(x, "x")
  bw <- attr(x, "bw") 
  
  inf <- min(xx)
  sup <- max(xx)
  z <- seq(inf, sup, (sup - inf)/1024)
  lz <- x(z)

  if (is.null(ylim)) ylim <- range(lz)    
  if (is.null(main)) main <- "Empirical Lientz's function"
  if (is.null(xlab)) xlab <- "x"
  if (is.null(ylab)) ylab <- "Sn(x)"
    
  plot.default(z, lz, main = main, xlab = xlab, ylab = ylab, ylim = ylim, ...)
  points.default(xx, rep(ylim[1],length(xx)), pch = "'", col = 4)
  legend("topleft",legend = c("Regular grid", "x"), col = c(1,4), pch = 19, bg = "white")
  
  if (zoom) {
    cat("you can zoom on the graph (press 'Esc' to escape)\n")  
    lc <- locator(2)
    while (!is.null(lc)) {
      xlim <- sort(c(lc$x[1], lc$x[2]))
      ylim <- sort(c(lc$y[1], lc$y[2]))
      plot.lientz(x, zoom = FALSE, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
      lc <- locator(2)
    }  
  }
  
  return(invisible(NULL))
}


print.lientz <-
function(x,             # an object of class 'lientz'
         digits = NULL,
         ...)
{
  if (!inherits(x, "lientz")) {
    stop("argument 'x' must inherit from class 'lientz'")
  }
  bw <- attr(x, "bw")
  call <- attr(x, "call")
  #xx <- attr(x, "x")
  cat("Empirical Lientz function\n")
  cat("Call:",deparse(call),"\n")
  cat("bw =", format(bw, digits = digits), "\n")
}


mlv.lientz <-
function(x,                       # sample (the data) or object of class 'lientz'
         bw = NULL,               # bandwidth
         abc = FALSE,            # if FALSE, 'optim' is used
         par = shorth(x),         # initial value used in 'optim'
         optim.method = "BFGS",   # method used in 'optim'
         ...)
{
#########################
# Lientz's mode estimator
#########################

  ## Initialization
  if (!inherits(x, "lientz")) {
    Sn <- lientz(x, bw)
  } else {
    Sn <- x
    x <- attr(Sn, "x")
  }
    
  if (!abc) {
    mini <- optim(par, fn = Sn, method = optim.method, control=list(fnscale=1),...)
    M <- mini$par
    attr(M, "value") <- mini$value
    attr(M, "counts") <- mini$counts
    attr(M, "convergence") <- mini$convergence
    attr(M, "message") <- mini$message
  } else {
    Sn <- Sn(x)
    M <- mean(x[Sn == min(Sn)])
  }
    
  ## Output
  return(M)
}
