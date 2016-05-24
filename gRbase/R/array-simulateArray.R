
## Simulate n observations from the array x conditional on
## the variables in margin (a vector of indices) takes values
## given by margin.value
##


simulateArray <- function(x, nsim=1, margin, value.margin){
  if(missing(margin)) {
    rhs       <- NULL
    lhs.dim   <- dim(x)
    lhs.names <- names(dimnames(x))
  } else {
    rhs       <- margin
    idx       <- (1:length(dim(x)))[-rhs]
    lhs.dim   <- (dim(x))[idx]
    lhs.names <- names(dimnames(x))[-rhs]
  }
  ##cat(sprintf("rhs=%s, lhs.dim=%s, lhs.names=%s\n",
  ##            toString(rhs), toString(lhs.dim), toString(lhs.names)))
  llhs.dim <- length(lhs.dim)
  pp   <- tableSlice(x, margin=rhs, level=value.margin)
  ##print(pp)
  samp <- sample(length(pp), size=nsim, replace=TRUE, prob=pp)
  ##print(samp)
  cp   <- cumprod(c(1, lhs.dim[-llhs.dim]))
  res  <- matrix(0, nrow=nsim, ncol=llhs.dim)
  for(j in 1:nsim){
    res[j,] <- 1+( (samp[j]-1) %/% cp ) %% lhs.dim
  }
  colnames(res) <- lhs.names
  res
}
