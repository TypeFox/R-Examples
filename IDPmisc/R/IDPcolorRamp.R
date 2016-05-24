### IDPcolorRamp.R

IDPcolorRamp <- function (n,
                          colInt = data.frame(
                            h = c(0.47, 0.28, 0.16, 0.00, 1.00, 0.8),
                            s = c(0.31, 0.55, 0.7, 0.8, 0.8, 1.00),
                            v = c(1, 1, 1, 1, 1, 0.4)),
                          fr     = c(0.27, 0.27, 0.27, 0))

  ## Author: Rene Locher
  ## Version: 2005-10-17
{
  if(!is.vector(n)) stop("n must be a single value\n")
  if(length(n)!=1) stop("n must be a single value\n")
  if(n<1) stop("n must be >= 1\n")
  n <- ceiling(n) ## n must be an integer
  nsr <- length(fr)+1
  if(sum(fr)<=0 | sum(fr)>1 | any(fr<0))
    stop("fr must sum to x with  0 < x < 1 where all fr are positive\n")
  if(nrow(colInt)!=nsr+1)
    stop("nrow(colInt) must be equal to length(fr)+2\n")
  if (nsr<2) stop("Number of Color Subramps must be at least 2\n")

  col <- data.frame(srn=1:nsr,ncol=rep(0,nsr),mod=rep(0,nsr))
  col$ncol[-nsr] <- trunc(fr*n)

  ## first subramp must have at least 1 color
  if(col$ncol[1]<1) col$ncol[1] <- 1

  col$mod <- c(fr,0)*n-col$ncol  ## ncol of last colRamp still 0
  delta <- trunc(sum(col$mod))

  if(delta>0) {                  ## distribute the remaining colors
    col1 <- col[-c(1,nsr),]
    ii <- order(col1$mod,decreasing = TRUE)
    ii <- ii[1:delta]
    col1$ncol[ii] <- col1$ncol[ii]+1
    col[-c(1,nsr),"ncol"] <- col1$ncol
    col$mod <- c(fr,0)*n-col$ncol
  }

  col$ncol[nsr] <- n-sum(col$ncol)

  delta <- col$ncol[nsr]-1
  if(delta<0) { ## ncol of last colRamp must be at least 1
    delta <- abs(delta)
    col1 <- col[-c(1,nsr),]
    ii <- order(col1$mod,col1$srn)
    ii <- ii[1:delta]
    col1$ncol[ii] <- col1$ncol[ii]-1
    col[-c(1,nsr),"ncol"] <- col1$ncol
    col$mod <- c(fr,col$ncol/n)*n-col$ncol
  }

  colRamp <- rep(NA,n)

  colRamp[1:col$ncol[1]] <-
    hsv(h = seq(colInt[1,1], colInt[2,1], length.out = col$ncol[1]),
        s = seq(colInt[1,2], colInt[2,2], length.out = col$ncol[1]),
        v = seq(colInt[1,3], colInt[2,3], length.out = col$ncol[1]))

  for (ii in 2:nsr) {
    if(col$ncol[ii]>0)
      colRamp[(sum(col$ncol[1:(ii-1)])+1):sum(col$ncol[1:ii])] <-
        hsv(h = seq(colInt[ii,1], colInt[ii+1,1],
              length.out = col$ncol[ii]+1),
            s = seq(colInt[ii,2], colInt[ii+1,2],
              length.out = col$ncol[ii]+1),
            v = seq(colInt[ii,3], colInt[ii+1,3],
              length.out = col$ncol[ii]+1))[-1]
  }
  return(colRamp)
} ## IDPcolorRamp

