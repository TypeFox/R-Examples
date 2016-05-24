subsets <- function(k) {
  subsets <- expand.grid(rep(list(0:1),k))
  ad  <- apply(subsets,1,sum)
  subsets[order(ad)[-1],]
}
