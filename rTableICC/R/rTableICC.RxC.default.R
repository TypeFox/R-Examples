rTableICC.RxC.default <-
function(p=NULL,theta,M,row.margins=NULL,col.margins=NULL,sampling="Multinomial",N=1,lambda=NULL,zero.clusters=FALSE,print.regular=TRUE,
         print.raw=FALSE){
  
  options(warn=-1)

  check(p, theta, M, N, K=NULL, row.margins, col.margins, lambda, sampling, ICC=TRUE, structure="RxC")

  tbl = rtableICC.RxC.main(p,theta,M,row.margins,col.margins,sampling,N,lambda,zero.clusters,print.regular,print.raw)
  tbl$call = match.call()
  class(tbl) = c("rTableICC.RxC", "rTableICC")
  tbl
}
