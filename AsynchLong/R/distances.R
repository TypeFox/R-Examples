distanceLV <- function(xIs, 
                       data.x, 
                       data.y, ...) {

  dis <- list()

  for( i in 1L:nrow(data.y) ) {

    dis[[i]] <- 0.0

    if( xIs[[i]]$n < 0.5 ) next

    xtime <- data.x[xIs[[i]]$v,2L]
    xtime[ xtime > data.y[i,2L]] <- NA

    tst <- which.max(xtime)

    if( length(tst) < 0.5 ) {
      xIs[[i]]$v <- array(dim=0L)
      xIs[[i]]$n <- 0L
      next
    }

    dis[[i]] <- 1.0

    xIs[[i]]$v <- xIs[[i]]$v[tst]
    xIs[[i]]$n <- 1L

  }

  return(list("dis" = dis, 
              "xIs" = xIs))

}


distanceTI <- function(xIs, 
                       data.x, 
                       data.y, 
                       bandwidth, 
                       kType, ...) {

  dis <- list()

  nCov <- ncol(data.x) - 2L

  lbw <- length(bandwidth)

  spr <- isTRUE(all.equal(lbw,nCov))

  for( i in 1L:nrow(data.y) ) {

    if( xIs[[i]]$n < 0.5 ) next

    mat <- matrix(data = 0.0, 
                  nrow = xIs[[i]]$n, 
                  ncol = nCov)

    kernelTime <- data.y[i,2L] - data.x[xIs[[i]]$v,2L]

    for( k in 1:lbw ) {
      mat[,k] <- local_kernel(t = kernelTime, 
                              h = bandwidth[k], 
                              kType = kType)
    }
    if( !spr ) mat[,2L:nCov] <- mat[,1L]

    dis[[i]] <- mat

  }

  return(list("dis" = dis, 
              "xIs" = xIs))

}


distanceTD <- function(xIs, 
                       data.x, 
                       data.y, 
                       bandwidth, 
                       kType,
                       tt, ...) {

  lbw <- length(bandwidth)

  if( !isTRUE(all.equal(lbw,(ncol(data.x)-2L))) ) {
    bandwidth <- rep(bandwidth, ncol(data.x)-2L)
  }

  lbw <- length(bandwidth)

  funcXY <- function(x, dfx, dfy) {

    kernelTime <- c(tt - dfx[,2L],tt - dfy[,2L])

    return(local_kernel(t = kernelTime, 
                        h = bandwidth[x], 
                        kType = kType))

  }

  nx <- nrow(data.x)
  ny <- nrow(data.y)
  
  disXY <- apply(X = array(1L:lbw), 
                 MARGIN = 1L, 
                 FUN = funcXY,  
                 dfx = data.x,  
                 dfy = data.y)
  
  disX <- disXY[1:nx,,drop=FALSE]
  disY <- disXY[(nx+1):(nx+ny),,drop=FALSE]
  
  dis <- lapply(X = 1L:ny, 
                FUN = function(x){
                        disY[x,]*disX[xIs[[x]]$v,,drop=FALSE]
                      })
  
  return(list("dis" = dis, 
              "xIs" = xIs))

}

distanceWLV <- function(xIs, 
                        data.x, 
                        data.y,
                        bandwidth, 
                        kType, ...) {

  dis <- list()

  nCov <- ncol(data.x) - 2L

  lbw <- length(bandwidth)

  spr <- isTRUE(all.equal(lbw,nCov))

  for( i in 1L:nrow(data.y) ) {

    dis[[i]] <- 0.0

    if( xIs[[i]]$n < 0.5 ) next

    xtime <- data.x[xIs[[i]]$v,2L]
    xtime[ xtime > data.y[i,2L]] <- NA

    tst <- which.max(xtime)
    if( length(tst) < 0.5 ) {
      xIs[[i]]$v <- array(dim=0L)
      xIs[[i]]$n <- 0L
      next
    }

    xtime <- xtime[tst]

    xIs[[i]]$v <- xIs[[i]]$v[tst]
    xIs[[i]]$n <- 1L

    kernelTime <- data.y[i,2L] - xtime

    mat <- matrix(data = 0.0, 
                  nrow = 1L, 
                  ncol = nCov)

    for( k in 1:lbw ) {
      mat[,k] <- local_kernel(t = kernelTime, 
                              h = bandwidth[k], 
                              kType = kType)
    }
    if( !spr ) mat[,2L:nCov] <- mat[,1L]

    dis[[i]] <- mat

  }

  return(list("dis" = dis, 
              "xIs" = xIs))

}


distanceHK <- function(xIs, 
                       data.x, 
                       data.y, 
                       bandwidth, 
                       kType, ...) {

  dis <- list()

  nCov <- ncol(data.x) - 2L

  lbw <- length(bandwidth)

  spr <- isTRUE(all.equal(lbw,nCov))

  for( i in 1L:nrow(data.y) ) {

    if( xIs[[i]]$n < 0.5 ) next

    xtime <- data.x[xIs[[i]]$v,2L]
    tst <- xtime < data.y[i,2L]

    mat <- matrix(data = 0.0, 
                  nrow = xIs[[i]]$n, 
                  ncol = nCov)

    kernelTime <- data.y[i,2L] - data.x[xIs[[i]]$v,2L]

    for( k in 1L:lbw ) {
      mat[,k] <- local_kernel(t = kernelTime, 
                              h = bandwidth[k], 
                              kType = kType)
    }
    if( !spr ) mat[,2L:nCov] <- mat[,1L]

    dis[[i]] <- mat*tst

  }

  return(list("dis" = dis, 
              "xIs" = xIs))

}

