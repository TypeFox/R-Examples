`fitweibull4` <-
function(x, y = NULL, p0 = c(0.1, 50, 5, 100), linint = -1, maxit=1000) {
  ## local functions: optimization criteria
  fmin <- function(p, x, y) sum((fweibull4(x, p) - y)^2)
  mse  <- function(x, y)    mean((x - y)^2)
  
  ## the following allows that both columns may be in x
  xy <- xy.coords(x, y)
  x  <- xy$x
  y  <- xy$y

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
    

  ## Heuristics to find start parameters
  ## ... not implemented for Weibull 4
  ## please use appropriate p0 vector

  ## optimize
  wrn <- options(warn = -1) # switch warnings off
  tno <- optim(p0, fmin, method="BFGS",
               # method="L-BFGS-B",
               #lower = c(0, 0, 0, 0),
               #upper = c(1, 1000, 1000, 1000),
               x = xx, y = yy,
               control = list(maxit = maxit))
  options(wrn) ## return to default options
  p <- as.numeric(tno$par)
  r2 <- 1 - var(fweibull4(x, p) - y)/var(y)
  ret <- list(p    = p,
       ymax = max(y),
       r2   = 1 - var(fweibull4(x, p) - yrel)/var(yrel),
       fit  = data.frame(x = xx, y = yy, f = fweibull4(xx, p))
  )
  class(ret) <- c("list", "cardiFit", "cardi4Fit")
  ret
}

