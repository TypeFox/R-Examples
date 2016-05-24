wccmap <- function(x, obj) {
  codes <- x$codes
  trwdth <- x$trwdth
  wccs <- rep(0, nrow(codes))
  if (trwdth > 0)
    wghts <- 1 - (0:trwdth)/trwdth
  else
    wghts <- 1
  
  acor <- wac(obj, trwdth, wghts)
  for (i in 1:nrow(codes))
    wccs[i] <- wcc(obj, codes[i,], trwdth, wghts,
                   c(acor, x$acors[i]))

  wccs
}

