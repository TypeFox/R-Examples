smoothColors<-function(...,alpha=NA){
 args <- list(...)
 r <- g <- b <- NULL
 while(length(args) > 0) {
  if(!is.character(args[[1]]))
   stop("Usage: smoothColors(\"color name\",[n|\"color name\"],...,\"color name\")")
  arglen<-length(args)
  if(arglen > 1){
   if(is.numeric(args[[2]])){
    lastarg<-2
    # args[[lastarg]] should be a color name
    while(is.numeric(args[[lastarg]])) {
     lastarg<-lastarg+1
     # make sure that there are enough arguments left
     if(lastarg > arglen) stop("bad argument list")
    } 
    ## do interpolate:
    from <- col2rgb(args[[1]])
    too <- col2rgb(args[[lastarg]])
    ## generate args[[2]] colors between specified colors:
    n <- args[[2]]+2 # add 2 for start and finish
    ## chop off last one since it will be added on the next iteration:
    r <- c(r,seq(from[1,],too[1,],length=n))
    i <- length(r)
    r <- r[-i]
    g <- c(g,seq(from[2,],too[2,],length=n))
    g <- g[-i]
    b <- c(b,seq(from[3,],too[3,],length=n))
    b <- b[-i]
    ## cut color and n from list and back we go
    args <- args[-(1:(lastarg-1))]
   }
   else {
    ## insert color, chop off 1
    cc <- col2rgb(args[[1]])
    r <- c(r,cc[1,])
    g <- c(g,cc[2,])
    b <- c(b,cc[3,])
    args <- args[-1]
   }
  }
  else {
   ## insert color, chop off 1
   cc <- col2rgb(args[[1]])
   r <- c(r,cc[1,])
   g <- c(g,cc[2,])
   b <- c(b,cc[3,])
   args <- args[-1]
  }
 }
 if(is.na(alpha)) rgb(r,g,b,maxColorValue=255)
 else rgb(r,g,b,alpha=alpha,maxColorValue=255)
}
