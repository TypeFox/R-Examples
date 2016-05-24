#' @export 
print.IntIndex <- function(x,...){
  mlist <- split(x$Mindex,rep(1:length(x$Mstrata),diff(c(0,x$Mstrata))))
  p <- x$petoInt[1,]
  q <- x$petoInt[2,]
  pqnames <- paste("(p;q)=",paste("(",p,";",q,"]",sep=""))
  pqnames[p==q] <- paste("(p;q)=",paste("[",p[p==q],";",q[p==q],"]",sep=""))
  names(mlist) <- pqnames
  Mlist <- lapply(mlist,function(u){
    L <- x$obsInt[1,u]
    R <- x$obsInt[2,u]
    out <- paste("(",L,";",R,"]",sep="")
    out[L==R] <- paste("[",L[L==R],";",R[L==R],"]",sep="")
    out
  })
  print(Mlist)
  Ilist <- split(x$Iindex,rep(1:length(x$Istrata),diff(c(0,x$Istrata))))
}
