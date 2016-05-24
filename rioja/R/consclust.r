chclust <- function(d, method="coniss") {
  if (!("dist" %in% class(d)))
     stop("Input must be a distance matrix")
  x <- as.matrix(d)
  if (!is.numeric(d))
     stop("Input matrix must be numeric")
  if (any(is.na(d)))
     stop("Missing values in input data")
  METHODS <- c("conslink", "coniss")
  method <- pmatch(method, METHODS)
  if(is.na(method))
     stop("Invalid clustering method")
  if(method == -1)
     stop("Ambiguous clustering method")
  ret <- .Call("chclust", x, as.integer(method), PACKAGE="rioja")
  if (is.character(ret))
     stop(ret)
  merge <- .find.groups(ret)
  tree <- list(merge=merge[, 1:2], height=sort(ret), seqdist = merge[, 3], order=1:nrow(x), labels=attr(d, "Labels"), method=METHODS[method], call=match.call(), dist.method = attr(d, "method"))
  class(tree) <- c("chclust", "hclust")
  tree
}

.find.groups <- function(height) {
  nr = length(height)
  x <- height
  merge <- matrix(nrow=nr, ncol=3)
  rec <- vector(mode="numeric", length=nr+1) 
  rec[] <- NA
  nG = 1
  for (i in 1:nr) {
     n <- which.min(x)
     minx <- min(x, na.rm=TRUE)
     merge[i, 3] <- minx
     if (is.na(rec[n]) & is.na(rec[n+1])) {
        merge[i,1] = -n
        merge[i,2] = -(n+1)
        rec[n] = nG
        rec[n+1] = nG
     } else {
        if (is.na(rec[n]) & !is.na(rec[n+1])) {
           merge[i,1] = -n
           merge[i,2] = rec[n+1]
           rec[n] = nG
           rec[rec == rec[n+1]] = nG
        } else {
           if (!is.na(rec[n]) & is.na(rec[n+1])) {
              merge[i,1] = rec[n]
              merge[i,2] = -(n+1)
              rec[rec == rec[n]] = nG
              rec[n+1] = nG
           } else {
              merge[i,1] = rec[n]
              merge[i,2] = rec[n+1]
              rec[rec == rec[n]] = nG
              rec[rec == rec[n+1]] = nG
           }
        }
     }
     x[n] <- NA
     nG <- nG+1
  }
  merge
}

plot.chclust <- function (x, labels = NULL, hang = 0.1,
              axes = TRUE, xvar=1:(length(x$height)+1), xlim=NULL, ylim=NULL, 
              x.rev = FALSE, y.rev=FALSE, horiz=FALSE, ...)
{
    merge <- x$merge
    if (!is.matrix(merge) || ncol(merge) != 2)
	  stop("invalid dendrogram")
    ## merge should be integer but might not be after dump/restore.
    if (any(as.integer(merge) != merge))
        stop("'merge' component in dendrogram must be integer")
    storage.mode(merge) <- "integer"
    n <- nrow(merge)
    height <- as.double(x$height)
    labels <-
	  if(missing(labels) || is.null(labels)) {
	     if (is.null(x$labels))
		      paste(1L:(n+1))
	     else
		      as.character(x$labels)
	  } else {
	     if(is.logical(labels) && !labels) # FALSE
		      character(n+1)
	     else
		      as.character(labels)
	  }
    if (is.null(xlim))
      xlim <- range(xvar)
    if (x.rev)
      xlim <- rev(xlim)
    if (is.null(ylim))
       ylim <- range(height)
    if (y.rev)
      ylim <- rev(ylim)
    if (horiz) {
      plot(0, type="n", xlim=ylim, ylim=xlim, xlab="", ylab="", axes=FALSE, ...)
    }
    else {
      plot(0, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", axes=FALSE, ...)
    }
    us <- par("usr")
    range <- ifelse(horiz, abs(us[2]-us[1]), abs(us[4]-us[3]))
    if (hang >= 0)
      hang <- range * hang
    else 
       hang <- -abs(us[3])
    if (horiz)
      ax.range <- c(us[3], us[4])
    else
      ax.range <- c(us[1], us[2])
    ccex <- par("cex")
    offset <- strwidth("m", units="figure", cex=ccex) * range
    if (horiz)
      .drawdend(n, 0, 0, merge, height, xvar, hang, labels, offset, srt=0, adj=ifelse(y.rev, 0, 1), horiz=horiz, ax.range=ax.range, ...)
    else 
      .drawdend(n, 0, 0, merge, height, xvar, hang, labels, offset, srt=90, adj=ifelse(y.rev, 0, 1), horiz=horiz, ax.range=ax.range, ...)
    if(axes & horiz)
        axis(1, ...) # , at=pretty(range(height)), ...)
    if(axes & !horiz)
        axis(2, ...) # , at=pretty(range(height)), ...)
    invisible()
}

# drawdend(int node, double *x, double *y, SEXP dnd_llabels, pGEDevDesc dd)
.drawdend <- function(node, x, y, merge, dnd_hght, dnd_xpos, dnd_hang, dnd_llabels, dnd_offset, horiz, ax.range, ...) {
# Recursive function for 'hclust' dendrogram drawing:
# Do left + Do right + Do myself
# "do" : 1) label leafs (if there are) and __
#        2) find coordinates to draw the  | |
#        3) return (*x,*y) of "my anchor"
#
# Essentially a R version of the C code of function drawdend in the base R.
# Modified to plot vertical or horizontal diagrams by Steve Juggins

    xx <- vector("numeric", length=4)
    yy <- vector("numeric", length=4)
    y = dnd_hght[node]
#    /* left part  */
    k = merge[node, 1]
    if (k > 0) {
       ret <- .drawdend(k, x, y, merge, dnd_hght, dnd_xpos, dnd_hang, dnd_llabels, dnd_offset, horiz, ax.range, ...)
       xl <- ret[1]
       yl <- ret[2]
   }
    else {
      u <- par("usr")
    	xl = dnd_xpos[-k]
     	yl = ifelse(dnd_hang >= 0, y - dnd_hang, min(-dnd_hang, ifelse(horiz, u[1], u[3])))
      yp <- yl-dnd_offset
      if (horiz) {
         if (xl >= min(ax.range) & xl <= max(ax.range)) {
            text(yp, xl, dnd_llabels[-k], xpd=NA, ...)
         }
      }
      else {
         if (xl >= min(ax.range) & xl <= max(ax.range)) 
           text(xl, yp, dnd_llabels[-k], xpd=NA, ...)    
      }
    }
#    /* right part */
    k = merge[node, 2]
    if (k > 0) {
       ret <- .drawdend(k, x, y, merge, dnd_hght, dnd_xpos, dnd_hang, dnd_llabels, dnd_offset, horiz, ax.range, ...)
       xr <- ret[1]
       yr <- ret[2]
    }
    else {
        u <- par("usr")
      	xr = dnd_xpos[-k]
      	yr = ifelse(dnd_hang >= 0, y - dnd_hang, min(-dnd_hang, ifelse(horiz, u[1], u[3])))
      	yp <- yr-dnd_offset
        if (horiz) {
          if (xr >= min(ax.range) & xr <= max(ax.range))
            text(yp, xr, dnd_llabels[-k], xpd=NA, ...)
        }
        else {
          if (xr >= min(ax.range) & xr <= max(ax.range)) 
            text(xr, yp, dnd_llabels[-k], xpd=NA, ...)
        }
    }
    xx[1] = xl; yy[1] = yl;
    xx[2] = xl; yy[2] = y;
    xx[3] = xr; yy[3] = y;
    xx[4] = xr; yy[4] = yr;
    if (horiz)
      lines(yy, xx, xpd=NA, ...)
    else
      lines(xx, yy, xpd=NA, ...)
    x = 0.5 * (xl + xr);
    return(c(x,y))
}
