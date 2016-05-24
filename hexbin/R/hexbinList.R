hexList <- function(x,y=NULL,given=NULL,xbins=30,shape=1,
                    xbnds = NULL, ybnds = NULL,
                    xlab = NULL, ylab = NULL)
{
  xl <- if (!missing(x)) deparse(substitute(x))
  yl <- if (!missing(y)) deparse(substitute(y))
  xy <- xy.coords(x, y, xl, yl)
   if(length(given)!=length(xy$x) | is.null(given))
    stop("Given is is different length from x and y")
  if(is.factor(given))
    given <- as.character(given)
  clss <- unique(given)
  if(is.null(xbnds))
    xbnds <- range(xy$x)
  if(is.null(ybnds))
    ybnds <- range(xy$y)  
  hbins <- vector(mode = "list",length=length(clss))
  i <- 1
  for(g in clss){
    hbins[[i]] <- hexbin(xy$x[given==g],xy$y[given==g],
                         xbins=xbins,shape=shape,xbnds=xbnds,ybnds=ybnds)
    i <- i+1
  }
  mx <- max(unlist(lapply(hbins,function(h)max(h@count))))
  mn <- min(unlist(lapply(hbins,function(h)min(h@count))))           
  hl <- new("hexbinList",n=length(hbins),hbins=hbins, Xbnds=xbnds,
            Ybnds=ybnds, Xbins=integer(xbins), Shape=shape, Bnames=clss,
            CntBnds=c(mn,mx))
  hl
}


setClass("hexbinList",
         representation(n="integer", hbins="vector",
                        Xbnds="numeric", Ybnds="numeric",
                        Xbins="numeric", Shape="numeric",
                        Bnames="character", CntBnds="numeric")

         )


bnds.check <- function(binlst, xb = TRUE, yb = TRUE)
{
    xb <-
        if(xb) {
            b <- binlst[[1]]@xbnds
            all(unlist(lapply(binlst, function(x, bnd) all(x@xbnds == bnd), b)))
        } else TRUE
    yb <-
        if(yb) {
            b <- binlst[[1]]@ybnds
            all(unlist(lapply(binlst, function(y, bnd) all(y@ybnds == bnd), b)))
        } else TRUE
    xb & yb
}

xbins.check <- function(binlst)
{
    xb <- binlst[[1]]@xbins
    all(unlist(lapply(binlst, function(y, xbin)all(y@xbins == xbin), xb)))
}

shape.check <- function(binlst)
{
    xs <- binlst[[1]]@shape
    all(unlist(lapply(binlst, function(y, xsh)all(y@shape == xsh), xs)))
}

list2hexList <- function(binlst)
{
  if(length(binlst) < 2)
    stop(" need at least 2 hex bin objects")
  if(!all(unlist(lapply(binlst, is, "hexbin"))))
    stop("All Elements of list must be hexbin objects")
  if(!bnds.check(binlst))
    stop("All bin objects in list need the same xbnds and ybnds")
  if(!xbins.check(binlst))
    stop("All bin objects in list need the same number of bins")
  if(!shape.check(binlst))
    stop("All bin objects in list need the same shape parameter")
  mx <- max(unlist(lapply(binlst,function(h)max(h@count))))
  mn <- min(unlist(lapply(binlst,function(h)min(h@count))))           
  xbins <- binlst[[1]]@xbins
  xbnds <- binlst[[1]]@xbnds
  ybnds <- binlst[[1]]@ybnds
  shape <- binlst[[1]]@shape
  hl <- new("hexbinList",n=length(binlst),hbins=binlst, Xbnds=xbnds,
            Ybnds=ybnds, Xbins=xbins, Shape=shape,
            Bnames=names(binlst), CntBnds=c(mn,mx))
  hl
}

setAs("list","hexbinList",function(from)list2hexList(from))

#setMethod("[", "hexbinList", function(hbl,i,...)
#{
#    if( length(list(...)) > 0 )
#        stop("extra subscripts cannot be handled")
#    if(missing(i)) hbl
#    hbl@hbins[i]
#})

##setMethod("[[", "hexbinList", function(hbl)
##{

##})
