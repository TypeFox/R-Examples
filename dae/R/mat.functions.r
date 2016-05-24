"elements" <- function(x, subscripts)
#a function that returns, in a vector, the elements of x specified by the subscripts.
#x is an array
#subscripts is a two dimensional object interpreted as elements by dimensions
{ if (!(class(subscripts) == "matrix" | class(subscripts) == "data.frame"))
    stop("subscripts must be in a matrix or data.frame")
  n <- nrow(subscripts)
  nsub <- ncol(subscripts)
  if (nsub == 2)
  { sapply(1:n, function(i)x[subscripts[i,1],subscripts[i,2]])
  }
  else
    sapply(1:n, function(i)eval(parse(text=paste("x[", paste(subscripts[i,], collapse=","), "]"))))
}

"mat.I" <- function(order)
{ diag(rep(1, order))
}

"mat.J" <- function(order)
{ n <- order*order
  matrix(rep(1, n), nrow=order, ncol=order)
}

"mat.ar1" <- function(order, rho)
#function to form the correlation matrix of size order with an ar1 pattern for
#correlation parameter rho
{ n <- order*order
  ar1 <- matrix(rep(rho, n), nrow=order, ncol=order)
  power <- abs(outer(1:order, 1:order, "-"))
  ar1 <- ar1^power
  return(ar1)
}

"mat.exp" <- function(coordinates, rho) 
  { order <- length(coordinates)
    n <- order * order
    mat <- matrix(rep(rho, n), nrow = order, ncol = order)
    rownames(mat) <- colnames(mat) <- coordinates
    power  <- abs(outer(coordinates,coordinates,'-'))
    mat <- mat^power
    return(mat)
  }


"mat.banded" <- function(x, nrow, ncol)
  #Function to form a banded matrix from a vector of values
  #- band 1 is the diagonal, band 2 the first subdiagonal and so on
{ nband <- length(x)
  if (nband > min(nrow, ncol))
    stop("Have supplied values for more than ",min(nrow, ncol) ," bands")
  matrix <- matrix(0, nrow=nrow, ncol=ncol)
  for (i in 1:nband)
    matrix[row(matrix)==col(matrix)+i-1 | row(matrix)+i-1 == col(matrix)] <- x[i]
  return(matrix)
}


"mat.dirprod" <- function(A, B)
#function create the direct product of A and B
{ rA <- nrow(A); cA <- ncol(A)
  rB <- nrow(B); cB <- ncol(B)
  Aexp <- A[rep(1:rA, each=rB), rep(1:cA, each=cB)]
  Bexp <- eval(parse(text=paste("cbind(", paste(rep("B", cA), collapse=","), ")")))
  Bexp <- eval(parse(text=paste("rbind(", paste(rep("Bexp", rA), collapse=","), ")")))
  Aexp*Bexp
}


"mat.dirsum" <- function(matrices)
  #Function to form the direct sum of a list of matrices
{ if (!is.list(matrices))
    stop("Must supply a list for matrices")
  if (!all(unlist(lapply(matrices, is.matrix))))
    stop("All elements of the matrices list must be matrices")
  nr <- lapply(matrices, nrow)
  nc <- lapply(matrices, ncol)
  m <- sum(unlist(nr))
  n <- sum(unlist(nc))
  dsum <- matrix(0, nrow=m, ncol=n)
  r1 <- r2 <- c1 <- c2 <- 0
  for (i in 1:length(matrices))
  { r1 <- r2 + 1
    c1 <- c2 + 1
    r2 <- r2 + nr[[i]]
    c2 <- c2 + nc[[i]]
    dsum[r1:r2, c1:c2] <- matrices[[i]]
  }
  return(dsum)  
}
