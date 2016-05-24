polylo <-
  function (x, degree = 1, monomial = FALSE) 
{
  if (degree >= 4) 
    warning("This is not a smart polynomial routine. You may get numerical problems with a degree of 4 or more")
  x <- as.matrix(x)
  dn <- dimnames(x)
  dd <- dim(x)
  np <- dd[2]
  ad <- rep(1, ncol(x))
### Used to have a x=scale(x)
### That messed up predictions on new data
### So we remove it
###  x=scale(x)
  if (np == 1) 
    monomial <- TRUE
  if (degree > 1) {
    if (monomial) {
      ad <- seq(degree)
      px <- x
      cc <- sapply(split(paste(diag(np)), rep(seq(np), 
                                              rep(np, np))), paste, collapse = "")
      tx <- x
      for (i in 2:degree) {
        px <- px * tx
        x <- cbind(x, px)
        cc <- c(cc, sapply(split(paste(diag(np) * i), 
                                 rep(seq(np), rep(np, np))), paste, collapse = ""))
      }
    }
    else {
      matarray <- array(x, c(dd, degree))
      for (i in 2:degree) matarray[, , i] <- x^i
      matarray <- aperm(matarray, c(1, 3, 2))
      x <- matarray[, , np, drop=TRUE]
      ad0 <- seq(degree)
      ad <- ad0
      ncol.mat0 <- degree
      ncol.x <- degree
      d0 <- paste(ad0)
      cc <- d0
      for (ii in seq(np - 1, 1)) {
        index0 <- rep(seq(ncol.mat0), ncol.x)
        index <- rep(seq(ncol.x), rep(ncol.mat0, ncol.x))
        newad <- ad0[index0] + ad[index]
        retain <- newad <= degree
        mat0 <- matarray[, , ii, drop = TRUE]
        if (any(retain)) 
          newmat <- mat0[, index0[retain]] * 
            x[, index[retain]]
        else newmat <- NULL
       ddn <- paste(d0[index0[retain]], cc[index[retain]], 
                     sep = "")
        zeros <- paste(rep(0, nchar(cc[1])), collapse = "")
        cc <- paste(0, cc, sep = "")
        d00 <- paste(d0, zeros, sep = "")
        x <- cbind(mat0, x, newmat)
        cc <- c(d00, cc, ddn)
        ad <- c(ad0, ad, newad[retain])
        ncol.x <- length(ad)
      }
    }
    if (!is.null(dn)) 
      dn[[2]] <- cc
    else dn <- list(NULL, cc)
    dimnames(x) <- dn
  }
  attr(x, "degree") <- ad
  x
}
