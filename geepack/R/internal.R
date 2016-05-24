crossutri <- function(wave) {
  n <- length(wave)
  if (n == 1) return(NULL)
  ans <- rep(0, n*(n-1)/2)
  k <- 1
  for (i in 1:(n-1))
    for (j in (i+1):n) {
      ans[k] <- paste(wave[i], wave[j], sep=":")
      k <- k + 1
    }
  ans
}

genZcor <- function(clusz, waves, corstrv) {
  if (corstrv == 1) return (matrix(0,0,0))
  crs <- clusz * (clusz - 1) / 2
  if (corstrv == 2 || corstrv == 3) {
    ans <-  matrix(1, length(clusz), 1)
    ##ans <-  matrix(1, sum(crs), 1)
    colnames(ans) <- c("alpha")
  }
  else {
    id <- rep(1:length(clusz), clusz)
    z1 <- unlist(lapply(split(waves, id), crossutri))
    z2 <- unlist(crossutri(1:max(clusz)))
    z <- factor(z1,levels=unique.default(z2))
    ans <- model.matrix(~z - 1)
    znames <- paste("alpha", z2, sep = ".")
    colnames(ans) <- znames
  }
  ans
}


genZodds <- function(clusz, waves, corstrv, ncat) {
  if (corstrv == 1) return (matrix(0,0,0))
  crs <- clusz * (clusz - 1) / 2
  c2 <- ncat * ncat
  if (corstrv == 2 | corstrv == 3) {
    ans <- matrix(1, sum(crs) * c2, 1)
    colnames(ans) <- c("alpha")
  }
  else {
    id <- rep(1:length(clusz), clusz)
    z1 <- unlist(lapply(split(waves, id), crossutri))
    z2 <- unlist(crossutri(1:max(clusz)))
    z <- factor(z1,levels=unique.default(z2))
    z <- model.matrix(~z - 1)
    ind <- gl(sum(crs), c2)
    ans <- z[ind,]
    colnames(ans) <- paste("alpha", 1:dim(ans)[2], sep=".")
  }
  ans
}
