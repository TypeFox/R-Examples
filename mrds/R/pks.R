# Distance sampling Kolmogorov-Smirnov test
pks <- function(Dn,n){

  diff <- 1
  p <- exp(-2*n*Dn^2)
  i <- 1
  while(abs(diff)>.0000001){
    i <- i+1
    diff <- (-1)^(i-1)*exp(-2*n*(i*Dn)^2)
    p <- p+diff
  }

  return(2*p)
}
