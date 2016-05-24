
DissimDistance <- function(x, y, tx=NULL, ty=NULL) {
  
  # If both temporal indices are missing, equal sampling is assumed and    both series begin and end in the same timestamp.
  if (is.null(tx) && is.null(ty)) {
    tx <- seq(0, 1, length.out=length(x))
    ty <- seq(0, 1, length.out=length(y))
  }
  
  # If the temporal index of one of the series is missing, an equal starting and ending point is assumed and the series is sampled constantly in that interval.
  if (is.null(tx)) {
    tx <- seq(ty[1], ty[length(ty)], length.out=length(x))
  }
  if (is.null(ty)) {
    ty <- seq(tx[1], tx[length(tx)], length.out=length(y))
  }
  
  
  if (class(try(DissimInitialCheck(x, y, tx, ty))) == "try-error") {
    return(NA)
  } else {
    
  # If the temporal indices are different a global index is calculated 
  # taking into account both indexes.
  if (length(tx) != length(ty)) {
    ind <- indglobal(tx, ty)
  } else {
    # If the temporal indices are different a global index is calculated.
    if(any(tx != ty)) {
      ind <- indglobal(tx, ty)
    } else {
    # If the temporal indices are equal the global index is equal to them.
      ind <- tx
    }
  }

  # These arrays tell us in which interval each point is in the global index.
  xglobal <- indpos(tx, ind)
  yglobal <- indpos(ty, ind)

  # The linear approximation for each interval is created for series x.
  coefx <- xfun(x, tx)
  coefx1 <- coefx[[1]][xglobal[- length(xglobal)]]
  coefx2 <- coefx[[2]][xglobal[- length(xglobal)]]

  # The linear approximation for each interval is created for series y.
  coefy <- xfun(y, ty)
  coefy1 <- coefy[[1]][yglobal[- length(yglobal)]]
  coefy2 <- coefy[[2]][yglobal[- length(yglobal)]]

  # The coefficients of the trinomial are calculated.
  a <- (coefx1 - coefy1) ^ 2 
  b <- 2 * (coefx1 * coefx2 + coefy1 * coefy2 - coefx1 * coefy2 - coefx2 * coefy1) 
  c <- (coefx2 - coefy2) ^ 2

  # The distance is calculated.
  D <- c(1:(length(ind) - 1)) * 0
  # If a==0, then b==0 and then the distance is calculated as follows:
  case1 <- which(a == 0)
  D[case1] <- sqrt(c[case1]) * diff(ind)[case1]

  # If a>0, for univariate case, root is always 0  
  # root <- b[case2]^2 - 4*a[case2]*c[case2]

  # So the distance is calculated as: 
  
  case2 <- which(a > 0)
  ind1 <- ind[1:(length(ind) - 1)]
  ind2 <- ind[2:length(ind)]
  aa <- a[case2]
  bb <- b[case2]
  cc <- c[case2]
  indd1 <- ind1[case2]
  indd2 <- ind2[case2]
  aux1 <- (2 * aa * indd1 + bb) / (4 * aa) * sqrt(aa * indd1 ^ 2 + bb * indd1 + cc) - (bb ^ 2 - 4 * aa * cc) / (8 * aa * sqrt(aa))
  aux2 <- (2 * aa * indd2 + bb) / (4 * aa) * sqrt(aa * indd2 ^ 2 + bb * indd2 + cc) - (bb ^ 2 - 4 * aa * cc) / (8 * aa * sqrt(aa))
  D[case2] <- aux2 - aux1

  # The total distance is calculated as the sum of all the pieces.
  d <- sum(D)

  return(d)
  }
}


# This function checks for initial errors.
DissimInitialCheck <- function(x, y, tx, ty) {
  
  if (! is.numeric(x) | ! is.numeric(y)) {
    stop('The series must be numeric', call.=FALSE)
  }
  if (! is.vector(x) | ! is.vector(y)) {
    stop('The series must be univariate vectors', call.=FALSE)
  }
  if (length(x) <= 1 | length(y) <= 1) {
    stop('The series must have a more than one point', call.=FALSE)
  }
  if (any(is.na(x)) | any(is.na(y))) {
    stop('There are missing values in the series', call.=FALSE)
  }     
    if (tx[1] != ty[1] | tx[length(tx)] != ty[length(ty)]) {
      stop('Series must begin and end in the same timestamp', call.=FALSE)
    }
    if (any(tx < 0) | any(ty < 0)) {
      stop('The temporal indices must be positive', call.=FALSE)
    }
    if (any(diff(tx) <= 0) | any(diff(ty) <= 0)) {
      stop('The temporal indices must be strictly increasing', call.=FALSE)
    }
  
    if (length(tx) != length(x) | length(ty) != length(y)) {
      stop('The length of the time indice must be equal to the length of the series', call.=FALSE)
  }
}

# This function creates a global index by fusing two indices.
indglobal <- function(tx, ty) {
  ind <- as.numeric(unique(c(tx, ty)))
  ind <- ind[order(ind)]
  return(ind)
}

# This function gives the location of the time points in a global index in a smaller subindex that is contained in it.
indpos <- function(subindex, globalindex) {
  aux <- globalindex
  aux[globalindex %in% subindex] <- c(1:length(subindex))
  aux[! globalindex %in% subindex] <- NA
  aux <- na.locf(aux)
  return(aux)
}


# This function calculates the linear approximation of a function
# for a given time index.
  xfun <- function(x, tx) {
  x <- as.numeric(x)
  coef1 <- diff(x) / diff(tx)
  coef2 <- x[1:(length(x) - 1)] - coef1*tx[1:(length(tx) - 1)]
  return(list(coef1, coef2))
}

