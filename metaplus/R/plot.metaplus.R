plot.metaplus <- function(x,...,extrameta=NULL) {
    
  if (!inherits(x, "metaplus"))
    stop("Use only with 'metaplus' xs.\n")
  if (!is.null(extrameta))
      if (any(!sapply(extrameta,inherits,what="metaplus")))
        stop("Use only with 'metaplus' xs.\n")

  if (x$justfit) stop("Cannot use with objects fitted with justfit=TRUE")
  
  if (is.null(extrameta)) nsummaries <- 1
  else nsummaries <- length(extrameta)+1
  forest(x=x$yi,sei=x$sei,slab=x$slab,ylim=c(-0.5-nsummaries,length(x$yi)+3),...)
  sumxi <- x$results[1,1]
  sumci.lb <- x$results[1,2]
  sumci.ub <- x$results[1,3]
  mlab <- x$label
  if (!is.null(extrameta)) {
    for (iadd in 1:length(extrameta)) {
      thex <- extrameta[[iadd]]
      sumxi <- c(sumxi,thex$results[1,1])
      sumci.lb <- c(sumci.lb,thex$results[1,2])
      sumci.ub <- c(sumci.ub,thex$results[1,3])
      mlab <- c(mlab,thex$label)
     }
  }
  extravars = list(...)
  if (length(extravars)==0) addpoly(x=sumxi, ci.lb=sumci.lb,ci.ub=sumci.ub,mlab=mlab)
  else {
# only pass relevent parameters to addpoly
    extravars <- list(transf=extravars$transf,atransf=extravars$atransf,targs=extravars$targs,efac=extravars$efac,
                      cex=extravars$cex,digits=extravars$digits)
    extravars <- extravars[!sapply(extravars, is.null)]
    extravars <- c(list(x=sumxi, ci.lb=sumci.lb,ci.ub=sumci.ub,mlab=mlab),extravars)
    do.call("addpoly",extravars)
  }
  abline(h=0)
}