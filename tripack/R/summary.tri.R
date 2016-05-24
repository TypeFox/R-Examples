summary.tri<-function(object, ...)
{
  if(!inherits(object,"tri"))
    stop("object must be of class \"tri\"")
  tnabor<- integer(object$tlnew)
  nnabs <- integer(object$n)
  nptr <- integer(object$n)
  nptr1 <- integer(object$n)
  nbnos <- integer(object$n)
  ans<-.Fortran("troutq",
                 as.integer(object$nc),
                 as.integer(object$lc),
                 as.integer(object$n),
                 as.double(object$x),
                 as.double(object$y),
                 as.integer(object$tlist),
                 as.integer(object$tlptr),
                 as.integer(object$tlend),
                 as.integer(6),
                 nnabs=as.integer(nnabs),
                 nptr=as.integer(nptr),
                 nptr1=as.integer(nptr1),
                 tnabor=as.integer(tnabor),
                 nbnos=as.integer(nbnos),
                 na=as.integer(0),
                 nb=as.integer(0),
                 nt=as.integer(0),
                 PACKAGE = "tripack")
  ans<-list(n=object$n,
            na=ans$na,
            nb=ans$nb,
            nt=ans$nt,
            nc=object$nc,
            call=object$call)
  class(ans)<-"summary.tri"
  ans
}
