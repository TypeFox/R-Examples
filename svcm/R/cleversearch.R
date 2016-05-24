cleversearch <- function(fn, lower, upper, ngrid, startvalue, logscale=TRUE, clever=TRUE, verbose=FALSE) {

  ##construct grid
  ndims <- length(lower)
  grid <- NULL
  for (i in 1:ndims) {
    if (logscale) {
      grid <- cbind(grid, 10^seq(lower[i], upper[i], length=ngrid))
    } else {
      grid <- cbind(grid, seq(lower[i], upper[i], length=ngrid))
    }
  }
  
  fmin <- Inf
  fn1 <- function(pnew) fn(pnew)
  
  if (clever) {
    
    ##initialize
    if (missing(startvalue)) {
      if (ndims == 1) { ##start at the lowest possible parameter
        index <- 1 
      } else {          ##start in the middle of grid
        index <- floor(ngrid/2) * rep(1, ndims) 
      }
    } else {
      index <- NULL
      for (i in 1:ndims) {
        tmp <- max(which(order(c(startvalue[i], grid[,i])) == 1) - 1, 1)
        index <- c(index, tmp)
      }
    }    
    par <- rep(0, ndims)
    for (i in 1:ndims) {
      par[i] <- grid[index[i], i]
    }
    lookup <- array(NA, rep(ngrid, ndims))
    
    ##search
    move <- 1
    nstep <- 0
    while (move) {
      
      move <- 0
      for (i in 1:ndims) {
        
        lookupi <- index
        
        for (j in (index[i] - 1):(index[i] + 1)) {
          
          j <- max(min(j, ngrid), 1)
          lookupi[i] <- j
          if (is.na(lookup[t(lookupi)])) {
            pnew <- par
            pnew[i] <- grid[j, i]
            fnew <- fn1(pnew)
            lookup[t(lookupi)] <- fnew
            nstep <- nstep + 1
          } else {
            fnew <- lookup[t(lookupi)]
          }
          if (fnew < fmin) {            
            fmin <- fnew
            index[i] <- j
            par <- pnew
            move <- move + 1
            if (verbose == TRUE) {
              cat(paste("\nIndex: ", paste(index, collapse=","),
                        ", Moved in step: ", nstep, ", Objective: ", fmin, "\n",
                        sep=""))
            }            
          }
        }##j
      }##i
    }##while
    
    
  } else { ##full grid evaluation
    
    if (ndims==1) {
      nstep <- 0
      pnew <- rep(0, ndims)
      for (i in 1:ngrid) {
        pnew <- grid[i, 1]
        fnew <- fn1(pnew)
        nstep <- nstep + 1
        if (fnew < fmin) {
          fmin <- fnew
          par <- pnew
          if (verbose == TRUE) {
            cat(paste("\nIndex: ", i, ", Moved in step: ", nstep, ", Objective: ",
                      fmin, "\n", sep=""))
          }
        }
      }
    } else if (ndims==2){
      nstep <- 0
      pnew <- rep(0, ndims)
      for (i in 1:ngrid) {
        pnew[1] <- grid[i, 1]
        for (j in 1:ngrid) {
          pnew[2] <- grid[j, 2]    
          fnew <- fn1(pnew)
          nstep <- nstep + 1
          if (fnew < fmin) {
            fmin <- fnew
            par <- pnew
            if (verbose == TRUE) {
              cat(paste("\nIndex: ", paste(c(i,j), collapse=","),
                        ", Moved in step: ", nstep, ", Objective: ", fmin, "\n",
                        sep=""))
            }
          }
        }
      }      
    } else if (ndims==3){
      nstep <- 0
      pnew <- rep(0, ndims)
      for (i in 1:ngrid) {
        pnew[1] <- grid[i, 1]
        for (j in 1:ngrid) {
          pnew[2] <- grid[j, 2]
          for (k in 1:ngrid) {
            pnew[3] <- grid[k, 3]          
            fnew <- fn1(pnew)
            nstep <- nstep + 1
            if (fnew < fmin) {
              fmin <- fnew
              par <- pnew
              if (verbose == TRUE) {
                cat(paste("\nIndex: ", paste(c(i,j,k), collapse=","),
                          ", Moved in step: ", nstep, ", Objective: ", fmin, "\n",
                          sep=""))
              }
            }
          }
        }
      }
    } else {
      stop("full grid evaluation only available on 1d, 2d and 3d parameters!")
    }
  }
  return(list(par=par, value=fmin, counts=nstep))  
}

