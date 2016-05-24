multicomp.mmc <- function(x,
                          focus=dimnames(attr(x$terms,"factors"))[[2]][1],
                          comparisons="mca",
                          lmat,
                          lmat.rows=lmatRows(x, focus),
                          lmat.scale.abs2=TRUE,
                          ry,
                          plot=TRUE,
                          crit.point,
                          iso.name=TRUE,
                          estimate.sign=1,
                          x.offset=0,
                          order.contrasts=TRUE,
                          main,
                          main2,
                          focus.lmat,
                          ...) {
  ## if.R(r={
    multicomp.lm <- NA ## make R-2.6.0dev happy
    stop("multicomp.mmc works only in S-Plus.  Use mmc in R.")
##   },
##        s={})

##   {
##     ## Save a copy of the data.frame in frame=0 to put it where multicomp.lm
##     ## needs to find it when the example data is run through Splus CMD check.
##     my.data.name <- as.character(x$call$data)
##     if (length(my.data.name)==0)
##       stop("Please provide an lm.object calculated with an explicit 'data=my.data.frame' argument.")
##     undo.it <- (!is.na(match(my.data.name, objects(0))))
##     if (undo.it) old.contents <- get(my.data.name, frame=0)
## ##    assign(my.data.name, get(my.data.name), frame=0)
##     my.data <- try(get(my.data.name))
##     if (class(my.data)=="Error")
##       my.data <- try(get(my.data.name, frame=sys.parent()))
##     if (class(my.data)=="Error")
##       stop("Please send me an email with a reproducible situation that got you here. (rmh@temple.edu)")
##     assign(my.data.name, my.data, frame=0)
##   }

##   ## pairwise differences
##   if (missing(crit.point)) {
##     mc.mca <- multicomp.lm(x, focus, ..., comparisons=comparisons, plot=FALSE)
##     crit.point <- mc.mca$crit.point
##   }
##   else
##     mc.mca <- multicomp.lm(x, focus, ..., comparisons=comparisons, plot=FALSE,
##                            crit.point=crit.point)
##   oldClass(mc.mca) <-  c("multicomp.hh", "multicomp")
##   mc.mca$focus <- focus

##   ## group means
##   mc.none <- multicomp.lm(x, focus, ..., comparisons="none", plot=FALSE,
##                           crit.point=crit.point)
##   oldClass(mc.none) <-  c("multicomp.hh", "multicomp")
##   mc.none$focus <- focus
##   mc.none$method <- mc.mca$method
##   mc.none$height <- mc.none$table[,"estimate"] * 2
##   if (length(unlist(list(...)$adjust)) > 1) {
##     warning("\nPlease verify that these two equivalent names of means
## are in the same order.  If not, then change the order of positions
## in the lmat.rows argument to match the groups column.\n")
##     tmp <- cbind(groups=names(mc.none$table[,"estimate"]),
##                  lmat.rows=dimnames(mc.mca$lmat[lmat.rows,])[[1]])
##     dimnames(tmp)[[1]] <- seq(nrow(tmp))
##     print(tmp, quote=FALSE)
##   }
##     if (length(mc.none$table[,"estimate"]) !=
##         nrow(mc.mca$lmat[lmat.rows,, drop=FALSE]))
##       stop("Please specify lmat.rows with multicomp.mmc on a design with more than one factor.")
##   mc.mca$height <- (mc.none$table[,"estimate"] %*%
##                     abs(mc.mca$lmat[lmat.rows,]))[1,]
##   if (estimate.sign != 0) mc.mca <- multicomp.reverse(mc.mca, estimate.sign)

##   ## user-specified lmat or focus.lmat
##   if (!missing(focus.lmat) || !missing(lmat)) {
##     if (missing(lmat))
##       lmat <- lmatContrast(mc.none$lmat, focus.lmat)
##     if (lmat.scale.abs2)
##       lmat <- sweep(lmat, 2, apply(abs(lmat[lmat.rows, , drop=FALSE]), 2, sum)/2, "/")
##     mc.lmat <- multicomp.lm(x, focus, ..., comparisons="none", plot=FALSE,
##                             crit.point=crit.point,
##                             lmat=lmat)
##     oldClass(mc.lmat) <-  c("multicomp.hh", "multicomp")
##     if (!is.null(mc.lmat$message)) stop(mc.lmat$message)
##     mc.lmat$focus <- focus
##     mc.lmat$method <- mc.mca$method
##     mc.lmat$height <- (mc.none$table[,"estimate"] %*% abs(lmat[lmat.rows,]))[1,]
##     if (estimate.sign != 0) mc.lmat <- multicomp.reverse(mc.lmat, estimate.sign)
##   }

##   {
##     ## restore frame=0
##     if (undo.it) assign(my.data.name, old.contents, frame=0)
##     else remove(my.data.name, frame=0)
##   }

##   ## result
##   result <- list(mca=mc.mca, none=mc.none)
##   if (!missing(focus.lmat) || !missing(lmat)) result$lmat <- mc.lmat
##   if (!missing(ry)) result$ry <- ry
##   oldClass(result) <- c("mmc.multicomp", "list")

##   if (order.contrasts) {
##     result$none <- multicomp.order(result$none)
##     result$mca <- multicomp.order(result$mca)
##     if (!missing(focus.lmat) || !missing(lmat))
##       result$lmat <- multicomp.order(result$lmat)
##   }

##   if (!missing(main))  result$main  <- main
##   if (!missing(main2)) result$main2 <- main2

##   ## plot
##   if (plot) plot.mmc.multicomp(result, iso.name=iso.name, x.offset=x.offset)

##   return(result)
}

"[.mmc.multicomp" <- function(x, ..., drop = TRUE) {
 result <- NextMethod("[")
 oldClass(result) <- oldClass(x)
 result
}

## source("c:/HOME/rmh/HH-R.package/HH/R/multicomp.mmc.R")
