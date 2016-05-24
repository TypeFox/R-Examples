rTableICC.2x2xK.default <-
function(p,theta,M,sampling="Multinomial",N=0,lambda=NULL,zero.clusters=FALSE,print.regular=TRUE,print.raw=FALSE){

    options(warn=-1)

    check(p, theta, M, N,K=NULL, row.margins=NULL, col.margins=NULL, lambda, sampling, ICC=TRUE, structure="2x2xK")
    
    tbl = rtableICC.2x2xK.main(p,theta,M,sampling,N,lambda,zero.clusters,print.regular,print.raw)
    tbl$call = match.call()
    class(tbl) = c("rTableICC.2x2xK", "rTableICC")
    tbl
}
