herlang.shape.all <- function(phnum, lbound, ubound) {
  herlang.shape.ipart(1, numeric(phnum), 1, phnum, lbound, ubound, list())
}

herlang.shape.ipart <- function(pos, shape, mini, res, lb, ub, result) {
  if (mini > res) {
    return(result)
  } else {
    for (i in mini:res) {
      shape[pos] <- i
      result <- herlang.shape.ipart(pos+1, shape, i, res-i, lb, ub, result)
    }
    if (lb <= pos && pos <= ub) {
      result <- c(result, list(shape[1:pos]))
    }
    return(result)
  }
}

herlang.shape.increment <- function(shape, ubound, phsize) {
  ## append
  if (length(shape) == ubound || sum(shape) == phsize) {
    return(list())
  }
  result <- c(1, shape)
  retval <- list(result)
  
  ## add
  for (i in unique(shape)) {
    tmp <- shape
    l <- which(tmp == i)
    n <- length(l)
    tmp[l[n]] <- tmp[l[n]] + 1
    retval <- c(retval, list(tmp))
  }
  return(retval)
}

herlang.shape.decrement <- function(shape, lbound) {
  if (length(shape)==1 && shape[1]==1) {
    return(list())
  }
  retval <- list()
  ## subtract
  for (i in unique(shape)) {
    tmp <- shape
    l <- which(tmp == i)
    n <- length(l)
    tmp[l[1]] <- tmp[l[1]] - 1
    tmp <- tmp[tmp != 0]
    if (length(tmp) >= lbound)
      retval <- c(retval, list(tmp))
  }
  return(retval)
}

