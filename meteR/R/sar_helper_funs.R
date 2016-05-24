## figure out areas
#' @export
.findAreas <- function(spp, abund, row, col, x, y, Amin, A0) {
    if(is.null(spp) | is.null(abund)) { # no data
        if(missing(Amin)) { # use row and col as nrow and ncol
            if(length(row) == 1 & length(col) == 1) {
                ## for simplicity make sure nrow <= ncol
                nrow <- min(row, col)
                ncol <- max(row, col)
                Amin <- A0/(nrow*ncol)
            } else {
                stop('row and col should be scalers and interpreted as number of rows and cols when A0 provided')
            }
        } else { # use Amin to figure out nrow and ncol
            nrow <- ncol <- floor(sqrt(A0/Amin))
            Amin <- A0/(nrow*ncol)
        }
        row <- col <- NULL
    } else { # data provided
        if(!is.null(x) & !is.null(y)) { # case where x,y data provided turn point data into row col data
            ## calculate areas, ranges, etc
            xrng <- diff(range(x))
            yrng <- diff(range(y))
            
            ## make sure row (=y) <= col (=x)
            if(xrng < yrng) {
                temp <- y
                y <- x
                x <- temp
                rng <- diff(range(x))
                yrng <- diff(range(y))
            }
            
            A0 <- xrng*yrng
            
            if(!is.null(row) & !is.null(col)) { # if row and col given use them to get nrow and ncol
                if(length(row) == 1 & length(col) == 1) {
                    ## for simplicity make sure nrow <= ncol
                    nrow <- min(row, col)
                    ncol <- max(row, col)
                } else {
                    ## if row and col are vectors doesn't make sense with point data
                    stop('if using x,y location data, row and col should be scalers indicating number of desired rows and columns')
                }
            } else { # if row and col not given use max extent of x,y data and Amin to make grid
                ## try to get nrow and ncol such that min cell is as close to square as possible
                ncell <- floor(A0/Amin)
                nrow <- floor(yrng/sqrt(Amin))
                ncol <- round(ncell/nrow)
            }
            
            ## now we have nrow and ncol, use those to make grid
            Amin <- A0/(nrow*ncol) # needs to be done even if Amin given cause fitting rows and cols
                                   # in ranges might have changed area slightly
            rowEndPoints <- seq(min(y)-.Machine$double.eps, max(y)+.Machine$double.eps, length=nrow+1)[-1]
            colEndPoints <- seq(min(x)-.Machine$double.eps, max(x)+.Machine$double.eps, length=ncol+1)[-1]
            rowcol <- apply(cbind(x, y), 1, function(X) {
                rowPos <- X[2] - rowEndPoints
                colPos <- X[1] - colEndPoints
                rowPos[rowPos > 0] <- min(rowPos) - 1
                colPos[colPos > 0] <- min(rowPos) - 1
                r <- which.max(rowPos)
                c <- which.max(colPos)

                return(c(r, c))
            })
            
            row <- rowcol[1, ]
            col <- rowcol[2, ]
            
        } else if(length(row) != length(spp) | length(col) != length(spp)) {
            stop('either row and column must be given for each spp or individual, or x,y coordinates given')
        } else { # case where row and col provided for each record
            ## make sure row and col IDs are 1:nrow and 1:ncol
            if(!is.numeric(row) | max(row) != length(unique(row))) row <- as.numeric(as.factor(row))
            if(!is.numeric(col) | max(col) != length(unique(col))) col <- as.numeric(as.factor(col))
            
            ## number of cells in each direction, transpose data if needed so always more columns
            if(max(row) > max(col)) {
                temp <- col
                col <- row
                row <- temp
            }
            
            ## get nrow and ncol
            nrow <- max(row)
            ncol <- max(col)
        }
    }

    ## figure out vector of sizes in units of cells; right now only doublings supported
    maxDoubling <- .calcMaxDoubling(floor(log(nrow*ncol) / log(2)), nrow)
    
    return(list(areas = 2^(0:maxDoubling), row=row, col=col, nrow=nrow, ncol=ncol, Amin=Amin, A0=A0))
}



#========================================================================
## funciton uses recursion to calculate maximum possible doubling for an area with 
## given minimum dimension
.calcMaxDoubling <- function(doub, mindim) {
  if(doub %% mindim != 0) doud <- doub - 1
  
  if(2^doub > 2*mindim^2) {
    doub <- .calcMaxDoubling(floor(log(2*mindim^2) / log(2)), mindim)
  }
  
  return(doub)
}

#========================================================================
## function to find neighbors and return groups of neighbors of given size `a'
.getNeighbors <- function(a, nr, nc) {
  foo <- function() {
    addToCol <- 1:(nc - c2 + 1)
    addToRow <- 1:(nr - c1 + 1)
    
    groups <- vector('list', length(addToRow) * length(addToCol))
    for(i in 1:length(addToRow)) {
      for(j in 1:length(addToCol)) {
        temp <- expand.grid(col=addToCol[j] + 0:(c2-1),
                            row=addToRow[i] + 0:(c1 - 1))
        groups[[i + length(addToRow)*(j-1)]] <- paste(temp[, 2], temp[, 1], sep=',')
      }
    }
    
    return(data.frame(group=rep(1:(length(addToRow) * length(addToCol)), 
                                each=a),
                      cells=unlist(groups), stringsAsFactors=FALSE))
  }
  
  c1 <- max((1:nr)[a %% (1:nr) == 0])
  c2 <- a/c1
  groups1 <- foo()
  
  if(c2 < nr & c1 != c2) {
    c1 <- c2
    c2 <- a/c1
    groups2 <- foo()
    groups2[, 1] <- groups2[, 1] + max(groups1[, 1])
    return(rbind(groups1, groups2))
  } else {
    return(groups1)
  }
}

#========================================================================
## function to get number of species in groups of cells
.getSppInGroups <- function(spp, abund, row, col, groups, endemics=FALSE) {
    cellID <- paste(row, col, sep=',')
    cellGroup <- groups$group[match(cellID, groups$cells)]
    # browser()
    sppByGroup <- tapply(abund, list(cellGroup, spp), sum)
    sppByGroup[is.na(sppByGroup)] <- 0
    sppByGroup[sppByGroup > 0] <- 1
    
    if(endemics) {
    	theseEndem <- colSums(sppByGroup) == 1
    	return(rowSums(sppByGroup[, theseEndem, drop=FALSE]))
    } else {
    	return(rowSums(sppByGroup))
    }
}

#========================================================================
## makeSSF is not vectorized and too slow so make special
## function to extract Pi(0) and make it vectorized over n0
.getPi0 <- function(n0,A,A0) {
  if(A/A0 == 0.5) {
    pi0 <- 1/(1+n0)
  } else if(A == A0) {
    pi0 <- 0
  } else {
    eq52 <- .useEq52(n0,A,A0)
    pi0 <- numeric(length(eq52))
    
    if(any(eq52)) {
      pi0[eq52] <- 1 - n0[eq52]/(n0[eq52] + A0/A)
    }
    
    if(any(!eq52)) {
      pi0[!eq52] <- sapply(n0[!eq52], function(n) {
        metePi(0, meteSSF(n0=n, A=A, A0=A0)$La, n0=n)
      })
    }
  }
  
  return(pi0)
}

#========================================================================
## makeSSF is not vectorized and too slow so make special
## function to extract Pi(n0) make it vectorized over n0
.getPin0 <- function(n0,A,A0) {
  if(A/A0 == 0.5) {
    pin0 <- 1/(1+n0)
  } else if(A == A0) {
    pin0 <- 1
  } else {
    eq52 <- .useEq52(n0,A,A0)
    pin0 <- numeric(length(eq52))
    
    if(any(eq52)) {
      pin0[eq52] <- ((n0[eq52]*A / (A0 + n0[eq52]*A))^n0[eq52]) / ((A0 + n0[eq52] * A) / A0)
    }
    
    if(any(!eq52)) {
      pin0[!eq52] <- sapply(n0[!eq52], function(n) {
        metePi(n, meteSSF(n0=n, A=A, A0=A0)$La, n0=n)
      })
    }
  }
  
  return(pin0)
}


#========================================================================
## FOR UPSCALING:
## the following two equations form a system that has to be solved to get upscaled
## species richness at an area Aup == 2*A0

## eq (8) from Harte et al. 2009 Ecol Lett and eq (7.70) from Harte 2011 Oxford Press
## solved for Sup so it can be plugged into eq (9)
.Sup <- function(beta, S0, N0) {
  exp(-beta) * (S0 + 2*N0*(1-exp(-beta)) / (exp(-beta) - exp(-beta*(2*N0+1))) * (1 - exp(-beta*2*N0) / (2*N0-1)))
}

## eq (9) from Harte et al. 2009 Ecol Lett and eq (7.71) from Harte 2011 Oxford Press
.eq9 <- function(beta, S0, N0) {
  ns <- 1:(2*N0)
  sapply(beta, function(b) {
    sum(exp(-b*ns) * (.Sup(b, S0, N0)/(2*N0) - 1/ns))
  })
}

#========================================================================
## figure out range in which to find solution for beta
.solRng <- function(S0, N0) {
	## empirically derived
	a0 <- -1
	a1 <- -1.38
	
	lwr <- exp(3*a0) * (N0/S0)^(1.5*a1)
	upr <- exp(0.1*a0)*2.5 * (N0/S0)^(0.9*a1)
	
	return(c(lwr=lwr, upr=upr))
 }

## find the root to the upscale constraint and return solution
#' @importFrom stats uniroot
.solveUpscale <- function(S0, N0) {
  beta <- uniroot(.eq9, .solRng(S0, N0), S0=S0, N0=N0, tol=.Machine$double.eps)$root
  Sup <- .Sup(beta, S0, N0)
  
  return(Sup)
}
