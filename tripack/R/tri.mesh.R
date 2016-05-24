tri.mesh <- function(x,y=NULL,duplicate="error")
{
  if(is.null(x))
     stop("argument x missing.")
  if(is.null(y)){
    x1<-x$x
    y1<-x$y
    if (is.null(x1) || is.null(y1))
      stop("argument y missing and x contains no $x or $y component.")
  }
  else{
    x1<-x
    y1<-y
  }

  n <- length(x1)
  if(length(y1)!=n)
    stop("length of x and y differ.")
  # handle duplicate points:
  xy <- paste(x1, y1, sep =",")
  i <- match(xy, xy)
  if(duplicate!="error")
    {
      if(duplicate!="remove" & duplicate!="error" & duplicate!="strip"){
        stop("possible values for \'duplicate\' are \"error\", \"strip\" and \"remove\"") 
      }
      else{
        if(duplicate=="remove")
          ord <- !duplicated(xy)
        if(duplicate=="strip")
          ord <- (hist(i,plot=FALSE,freq=TRUE,breaks=seq(0.5,max(i)+0.5,1))$counts==1)
        x1 <- x1[ord]
        y1 <- y1[ord]
        n <- length(x1)
      }
    }
  else
    if(any(duplicated(xy)))
      stop("duplicate data points")

  ans<-.Fortran("trmesh",
                as.integer(n),
                x=as.double(x1),
                y=as.double(y1),
                tlist=integer(6*n-12),
                tlptr=integer(6*n-12),
                tlend=integer(n),
                tlnew=as.integer(0),
                tnear=integer(n),
                tnext=integer(n),
                tdist=double(n),
                ier=as.integer(0),
                PACKAGE = "tripack")
  if(ans$ier==0)
    {
      tri.obj<-list(n=n,x=x1,y=y1,tlist=ans$tlist,tlptr=ans$tlptr,
                    tlend=ans$tlend,tlnew=ans$tlnew,
                    nc=0,lc=0,call=match.call())
  }
  else
    stop("error in trmesh")
                  
  class(tri.obj)<-"tri"
  invisible(tri.obj)
}
