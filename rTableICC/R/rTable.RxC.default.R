rTable.RxC.default <-
function(p,row.margins=NULL,col.margins=NULL,sampling="Multinomial",N,lambda=NULL,print.raw=FALSE){

  options(warn=-1)
  
  check(p, theta=NULL, M=NULL, N, K=NULL, row.margins, col.margins, lambda, sampling, ICC=FALSE, structure="RxC")
  
  tbl = rtable.RxC.main(p,row.margins,col.margins,sampling,N,lambda,print.raw)
  tbl$call = match.call()
  class(tbl) = c("rTable.RxC","rTableICC")
  tbl
}
