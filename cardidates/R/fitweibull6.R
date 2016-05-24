`fitweibull6` <-
function(x, y = NULL, p0 = NULL, linint = -1, maxit=2000) {
  ## local functions: optimization criteria
  fmin <- function(p, x, y) sum((fweibull6(x, p) - y)^2)
  mse  <- function(x, y)    mean((x - y)^2)

  ## the following allows that both columns may be in x
  xy <- xy.coords(x, y)
  x <- xy$x
  y <- xy$y

  xm   <- length(x)
  yrel <- y / max(y)
  
  if (linint < 0) { # apply interpolation heuristics
    ## extremely steep peak?
    fs   <- ifelse(length(which(yrel > 0.5)) < 5, 0.5, 1)
    ## sam = mean sampling rate
    sam  <- 2 * median(diff(x))
    ## if x distance < 10 & not steep & > 35 values
    if (sam < 10 & fs == 1 & xm > 35) {
       xx <- x
       yy <- yrel
    } else {
      xx <- seq(x[1], max(x))
      yy <- approx(x, yrel, xx)$y
    }
  } else if (linint == 0) { # no interpolation
      xx <- x
      yy <- yrel
  } else { # interpolation with given time step
      xx <- seq(x[1], max(x), linint)
      yy <- approx(x, yrel, xx)$y
  }
    
  isnull    <-  is.null(p0)
  tosmall   <- (length(p0) < 6)
  if (isnull | tosmall) {
    # print("heuristics for p0 used")
    ## Find start parameters
    x1   <- which.max(yrel[1:xm])  ### <--- xm obsolet??
    xa   <- x[x1] - 30
    if (xa < 50) xa <- x[x1] - 20
    xb   <- min(x[length(yrel)], x[x1] + 30)
    steep <- 5
    test <-  mse(yrel, fweibull6(x, c(1, xa, steep, 0, xb, steep)))
    steep <- c(10, 5, 1)[1 + (test >= 0.5) + (test > 1.5)]
    p0   <- c(1, xa, steep, 0, xb, steep)
  }
  xtno <- list(
    ## first optimizer: stochastic SANN
    tno1 = optim(p0, fmin, x = xx, y = yy, method="SANN", 
                 control=list(maxit=maxit)),
    ## second optimizer, constrained
    tno2 = optim(p0, fmin, method="L-BFGS-B",
                 lower = c(0,     1,   1e-6, 0,  1,   1e-6),
                 upper = c(1.5, 365, 100,    1,  365, 100), 
                 x = xx, y = yy,
                 control = list(maxit = maxit)),
    ## last resort, fixed heuristic start values
    tno3 = optim(c(1, 50, 50, 0.1, 150, 50), fmin, x = xx, y = yy, method="SANN", 
                 control=list(maxit=maxit)
    )
  )       
  best <- which.min(lapply(xtno, function(l) l$value))
  tno <- xtno[[best]]
                      
  #print(tno$par)
  ## improve with third optimizer
  p0 <- tno$par
  tno <- optim(p0, fmin, x = xx, y = yy, method="BFGS",
     control=list(maxit=maxit))
  # print(tno$par)       
  p <- as.numeric(tno$par)
  r2 <- 1 - var(fweibull6(x, p) - yrel)/var(yrel)

  ret <- list(p    = p,
       ymax  = max(y),
       convergence = tno$convergence,
       r2   = r2, 
       fit  = data.frame(x = xx, y = yy, f = fweibull6(xx, p))
  )
  class(ret) <- c("list", "cardiFit", "cardi6Fit")
  ret
}

