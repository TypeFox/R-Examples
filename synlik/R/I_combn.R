.combn <- function(x, m, fun = NULL, simplify = TRUE, ...)
{
  #       DATE WRITTEN: 14 April 1994          LAST REVISED:  10 July 1995
  #       AUTHOR:  Scott Chasalow
  #
  #       DESCRIPTION:
  #             Generate all combinations of the elements of x taken m at a time. 
  #             If x is a positive integer,  returns all combinations
  #             of the elements of seq(x) taken m at a time.
  #             If argument "fun" is not null,  applies a function given
  #             by the argument to each point.  If simplify is FALSE,  returns 
  #             a list; else returns a vector or an array.  "..." are passed 
  #             unchanged to function given by argument fun,  if any.
  #       REFERENCE:
  #             Nijenhuis, A. and Wilf, H.S. (1978) Combinatorial Algorithms for 
  #             Computers and Calculators.  NY:  Academic Press.
  #       EXAMPLES:
  #             > combn(letters[1:4], 2)
  #             > combn(10, 5, min)  # minimum value in each combination
  #             Different way of encoding points:
  #             > combn(c(1,1,1,1,2,2,2,3,3,4), 3, tabulate, nbins = 4)
  #             Compute support points and (scaled) probabilities for a
  #             Multivariate-Hypergeometric(n = 3, N = c(4,3,2,1)) p.f.:
  #             > table.mat(t(combn(c(1,1,1,1,2,2,2,3,3,4), 3, tabulate,nbins=4)))
  #
  if(length(m) > 1) {
    warning(paste("Argument m has", length(m), 
                  "elements: only the first used"))
    m <- m[1]
  }
  if(m < 0)
    stop("m < 0")
  if(m == 0)
    return(if(simplify) vector(mode(x), 0) else list())
  if(is.numeric(x) && length(x) == 1 && x > 0 && trunc(x) == x)
    x <- seq(x)
  n <- length(x)
  if(n < m)
    stop("n < m")
  e <- 0
  h <- m
  a <- 1:m
  nofun <- is.null(fun)
  count <- .nCm(n, m, 0.10000000000000002)
  out <- vector("list", count)
  out[[1]] <- if(nofun) x[a] else fun(x[a], ...)
  if(simplify) {
    dim.use <- NULL
    if(nofun) {
      if(count > 1)
        dim.use <- c(m, count)
    }
    else {
      out1 <- out[[1]]
      d <- dim(out1)
      if(count > 1) {
        if(length(d) > 1)
          dim.use <- c(d, count)
        else if(length(out1) > 1)
          dim.use <- c(length(out1), count)
      }
      else if(length(d) > 1)
        dim.use <- d
    }
  }
  i <- 2
  nmmp1 <- n - m + 1
  mp1 <- m + 1
  while(a[1] != nmmp1) {
    if(e < n - h) {
      h <- 1
      e <- a[m]
      j <- 1
    }
    else {
      h <- h + 1
      e <- a[mp1 - h]
      j <- 1:h
    }
    a[m - h + j] <- e + j
    out[[i]] <- if(nofun) x[a] else fun(x[a], ...)
    i <- i + 1
  }
  if(simplify) {
    if(is.null(dim.use))
      out <- unlist(out)
    else out <- array(unlist(out), dim.use)
  }
  out
}

