ELOS <- function(pt, tau) {
  if (!inherits(pt, "probtrans"))
    stop("'pt' must be a 'probtrans' object")
  if (pt$direction != "forward")
    stop("ELOS only works if direction==\"forward\" ")
  tt <- pt[[1]]$time
  if (missing(tau)) tau <- max(tt)
  if (tau<min(tt)) stop("tau should be larger than predt")
  tmat <- pt$trans
  K <- nrow(tmat)
  res <- matrix(NA, K, K)
  for (k in 1:K) {
    ptk <- pt[[k]]
    ptk <- subset(ptk, time<=tau)
    if (any(ptk[,(2:(K+1))]<0))
      warning("negative transition probabilities present; ELOS may not be well defined")
    ntk <- nrow(ptk)
    ptk <- rbind(ptk, ptk[ntk,])
    ptk$time[ntk+1] <- tau
    res[k,] <- apply(diff(ptk$time)*ptk[1:ntk,1+(1:K)], 2, sum)
  }
  rownames(res) <- paste("from", 1:K, sep="")
  colnames(res) <- paste("in", 1:K, sep="")
  return(res)
}
