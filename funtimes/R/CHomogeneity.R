CHomogeneity <- function(Bu, Bv, Alpha, Beta, Delta){
  M <- cbind(Bu,Bv)
  med <- apply(M, 2, median)
  MED <- matrix(med, length(med), length(med))
  tmp <- abs(MED - t(MED))/(Beta-Alpha) <= Delta
  all(tmp)
}