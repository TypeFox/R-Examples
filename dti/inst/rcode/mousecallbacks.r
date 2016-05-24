xprod <- function(a, b) 
    c(a[2]*b[3] - a[3]*b[2],
       a[3]*b[1] - a[1]*b[3],
       a[1]*b[2] - a[2]*b[1])

vlen <- function(a) sqrt(sum(a^2))

angle <- function(a,b) {
    dot <- sum(a*b)
    acos(dot/vlen(a)/vlen(b))
}
mouseTrackball <- function(button = 1, dev = rgl.cur() ) {
    width <- height <- rotBase <- NULL
    userMatrix <- list()
    cur <- rgl.cur()
    
    screenToVector <- function(x, y) {
      radius <- max(width, height)/2
      centre <- c(width, height)/2
      pt <- (c(x, y) - centre)/radius
      len <- vlen(pt)

      if (len > 1.e-6) pt <- pt/len

      maxlen <- sqrt(2)
      angle <- (maxlen - len)/maxlen*pi/2
      z <- sin(angle)
      len <- sqrt(1 - z^2)
      pt <- pt * len
      return (c(pt, z))
    }
    
    trackballBegin <- function(x, y) {
        vp <- par3d("viewport")
        width <<- vp[3]
        height <<- vp[4]
        cur <<- rgl.cur()
        for (i in dev) {
            if (inherits(try(rgl.set(i, TRUE)), "try-error")) dev <<- dev[dev != i]
            else userMatrix[[i]] <<- par3d("userMatrix")
        }
        rgl.set(cur, TRUE)
        rotBase <<- screenToVector(x, height - y)
    }
    
    trackballUpdate <- function(x,y) {
        rotCurrent <- screenToVector(x, height - y)
        angle <- angle(rotBase, rotCurrent)
        axis <- xprod(rotBase, rotCurrent)
        mouseMatrix <- rotationMatrix(angle, axis[1], axis[2], axis[3])
        for (i in dev) {
            if (inherits(try(rgl.set(i, TRUE)), "try-error")) dev <<- dev[dev != i]
            else par3d(userMatrix = mouseMatrix %*% userMatrix[[i]])
        }
        rgl.set(cur, TRUE)
    }
    
    for (i in dev) {
        rgl.set(i, TRUE)
        rgl.setMouseCallbacks(button, begin = trackballBegin, update = trackballUpdate, end = NULL)
    }
    rgl.set(cur, TRUE)
}
mouseInterp <- function(button = 1, dev = rgl.cur(), fn, init = 0, range = NULL, direction=c(1,0)) {
    cur <- rgl.cur()
    time <- init
    x0 <- width <- height <- NULL
    
    interpBegin <- function(x, y) {
    	vp <- par3d("viewport")
        width <<- vp[3]
        height <<- vp[4]
        x0 <<- sum(direction*c(x,y))
    }
        
    interpUpdate <- function(x,y) {
        time <<- init + (sum(direction*c(x,y)) - x0)/width
        if (!is.null(range)) time <<- clamp(time, range[1], range[2])
        for (i in dev) {
            if (inherits(try(rgl.set(i, TRUE)), "try-error")) dev <<- dev[dev != i]
            else par3d(fn(time))
        }
        rgl.set(cur, TRUE)
    }
    
    interpEnd <- function() {
        init <<- time
    }
    
    for (i in dev) {
        rgl.set(i, TRUE)
        rgl.setMouseCallbacks(button, begin = interpBegin, update = interpUpdate, end = interpEnd)
    }
    rgl.set(cur, TRUE)
}
mouseZoom <- function(button = 1, dev = rgl.cur()) 
    mouseInterp(button,dev=dev,fn=par3dinterp(times=c(-4,4)/4, zoom=c(10^(-4),10^4),method="linear"),
                      init=log10(par3d("zoom"))/4,range=c(-4,4)/4,direction=c(0,-1))
 
mouseFOV <- function(button = 1, dev = rgl.cur())
    mouseInterp(button,dev=dev,fn=par3dinterp(times=c(1,179)/180, FOV=c(1,179),method="linear"), 
                      init=par3d("FOV")/180, range = c(1,179)/180, direction=c(0,1))

clamp <- function(value, low, high)
{
  if (value < low) {
    warning( paste("value clamped to ",low) ); 
    result <- low
  }
  else if (value > high) {
    warning( paste("value clamped to ",high) );
    result <- high
  }
  else {
    result <- value
  }

  return (result);
}
