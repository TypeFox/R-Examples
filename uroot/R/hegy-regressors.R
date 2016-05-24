
hegy.regressors <- function(x)
{
  n <- length(x)
  S <- frequency(x)
  isSeven <- (S %% 2) == 0
  isSevenp2 <- 2 + isSeven

  ML <- x
  for (i in seq_len(S-1))
    ML <- cbind(ML, lag(x, -i))
  ML <- window(ML, end = end(x))
  ML <- sapply(seq.int(0, S-1), function(x, y, n) c(rep(NA,length.out=x), y[seq_len(n-x)]), y=x, n=n)

  ypi <- matrix(nrow = n, ncol = S)
  ypi[,1] <- rowSums(ML)
  if (isSeven)
    ypi[,2] <- ML %*% rep(c(-1, 1), len = S)
  seqS <- seq_len(S)
  j <- 0
  sinesign <- -1
  id <- seq.int(isSevenp2, S, 2)
  ref <- ceiling(S) / 4
  for (i in id)
  {
    j <- j + 1
    seqw <- seqS * (2 * pi * j / S)
    ypi[,i] <- ML %*% cos(seqw)
    ypi[,i+1] <- sinesign * ML %*% sin(seqw)
    if (j == ref)
      sinesign <- -1 * sinesign
  }

  ypi <- rbind(NA, ypi[-n,])
  colnames(ypi) <- paste("Ypi", seq_len(S), sep="")
  ypi[-seq_len(S),]
}
