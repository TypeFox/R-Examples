rTable.2x2xK.default <-
function(p=NULL,sampling="Multinomial",N=NULL,K=NULL,lambda=NULL,print.raw=FALSE){

  options(warn=-1)
  
  check(p, theta=NULL, M=NULL, N,K, row.margins=NULL, col.margins=NULL, lambda, sampling, ICC=FALSE, structure="2x2xK")
  
  tbl = rtable.2x2xK.main(p,sampling,N,K,lambda,print.raw)
  tbl$call = match.call()
  class(tbl) = c("rTable.2x2xK", "rTableICC")
  tbl
  
}
