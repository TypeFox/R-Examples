ck_compute <- function(n_mo, sk, p) {
  ck <- 2*log(sk + 2) + ifelse(sk > 0, sk*log(exp(1)*p/sk), 0)
  return(ck)
}


