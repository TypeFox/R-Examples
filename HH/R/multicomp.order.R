multicomp.order <- function(mca, sort.by="height", sort.order=NULL) {

  ## When (!is.null(sort.order)) ## use sort.order as the order
  ## When is.null(sort.order)    ## (the default), use sort.by
  
  ## When sort.by=="height", take the result of a multiple comparisons
  ## analysis from multicomp.mmc and sort the contrasts by the reverse
  ## order of the heights.  This provides a multicomp object that will
  ## be plotted by plot.multicomp in the same order that
  ## plot.multicomp.mmc uses.

  ## When sort.by=="estimate", sort the contrasts by the reverse order
  ## of the contrast estimates.  This provides the same order as the
  ## half-normal plot.

  if (is.null(sort.order)) {
    if (sort.by=="height") {
      if (is.null(mca$height)) return(mca)
      else sort.order <- rev(order(mca$height))
    }
    else if (sort.by=="estimate")
      sort.order <- rev(order(mca$table[,"estimate"]))
    else stop('"sort.by" must be one of c("height", "estimate").')
  }
  
  mca$table <- mca$table[sort.order,,drop=FALSE]
  mca$lmat <- mca$lmat[,sort.order,drop=FALSE]
  mca$height <- mca$height[sort.order]
  mca
}


## ## trace(multicomp.mmc, exit=browser)
## ## source(hh("splus.library/multicomp.mmc.s"))


multicomp.label.change <- function(x, old="adj", new="new", how.many=2)
  UseMethod("multicomp.label.change")

multicomp.label.change.multicomp <-
  function(x, old="adj", new="new", how.many=2) {

    old.new.length <- max(length(old), length(new))
    old      <- rep(old,      length=old.new.length)
    new      <- rep(new,      length=old.new.length)
    how.many <- rep(how.many, length=old.new.length)

    names <- dimnames(x$table)[[1]]

    for (j in seq(along=old))
      for (i in 1:how.many[j])
        names <- sub(old[j], new[j], names, fixed=TRUE)
 
    dimnames(x$table)[[1]] <- names
    dimnames(x$lmat)[[2]] <- names
    if (!is.null(x$height)) names(x$height) <- names
    x
  }

multicomp.label.change.mmc.multicomp <-
  function(x, old="adj", new="new", how.many=2) {
  for (i in c("mca","none"))
    x[[i]] <- multicomp.label.change(x[[i]], old, new, how.many=how.many)
  x
}
