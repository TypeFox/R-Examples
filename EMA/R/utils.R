myPalette <- function (low = "white", high = c("green", "red"), mid = NULL, k = 50){
    low <- col2rgb(low)/255
    high <- col2rgb(high)/255
    if (is.null(mid)) {
        r <- seq(low[1], high[1], len = k)
        g <- seq(low[2], high[2], len = k)
        b <- seq(low[3], high[3], len = k)
    }
    if (!is.null(mid)) {
        k2 <- round(k/2)
        mid <- col2rgb(mid)/255
        r <- c(seq(low[1], mid[1], len = k2), seq(mid[1], high[1], len = k2))
        g <- c(seq(low[2], mid[2], len = k2), seq(mid[2], high[2], len = k2))
        b <- c(seq(low[3], mid[3], len = k2), seq(mid[3], high[3], len = k2))
    }
    rgb(r, g, b)
}

as.colors <- function(x, col.na="#E6E6E6", palette="rainbow", ...){
  v2m=FALSE

  if (is.na(palette) || is.null(palette)){
    stop("Palette must be specified")
  }
  
  if (is.vector(x)){
    v2m=TRUE
    x <- matrix(x, nrow=1)
  }
 
  xfactor <- factor(as.vector(as.matrix(x)))
  col.sel <- do.call(palette,list(n=length(levels(xfactor)), ...))

  colorlab<-apply(x,1,function(xx){
    for (i in levels(factor(xx))){
      ind<-which(levels(xfactor)==i)
      xx[which(xx==i)] <- col.sel[ind]
    }
    return (xx)
  })

  
  colorlab[which(is.na(colorlab))] <- col.na
  
  if (v2m){
    return(as.vector(colorlab))
  }
  else if(is.matrix(x) && is.vector(colorlab)){
    colorlab<-as.matrix(colorlab,ncol=ncol(x), nrow=nrow(x))
  }
  else{
    return(t(colorlab))
  }
}


intersectg <- function(...) {
   args <- list(...)
   emptyin <- unlist(lapply(args,is.null))
   if (any(emptyin))
       args <- args[which(!emptyin)]
   nargs <- length(args) 
   if(nargs <= 1) {
     if(nargs == 1 && is.list(args[[1]])) {
       do.call("intersectg", args[[1]])
     } else {
       stop("cannot evaluate intersection fewer than 2 arguments")
     }
   } else if(nargs == 2) {
     intersect(args[[1]], args[[2]])
   } else {
     intersect(args[[1]], intersectg(args[-1]))
   }
}

setdiffg <- function(...) {
    args <- list(...)
    emptyin <- unlist(lapply(args,is.null))
    if (any(emptyin))
        args <- args[which(!emptyin)]
    nargs <- length(args)
    if(nargs <= 1) {
        if(nargs == 1 && is.list(args[[1]])) {
            do.call("setdiffg", args[[1]])
        } else {
            stop("cannot evaluate intersection fewer than 2 arguments")
        }
    } else if(nargs == 2) {
        setdiff(args[[1]], args[[2]])
    } else {
        setdiff(setdiffg(args[-length(args)]), args[[length(args)]])
    }
}
