#' Plots linkage maps
#' 
#' Plot linkage map (either as input object or as stored in mpcross object). 
#' Can also highlight QTL regions when used with qtlmap function. 
#' @export
#' @param object Either \code{mpcross} or \code{map} object
#' @param chr Chromosomes to plot
#' @param max.dist Plotting paramters. See \code{\link[wgaim]{link.map.cross}}
#' @param marker.names Whether to plot marker names
#' @param tick Plotting parameters. See \code{\link[wgaim]{link.map.cross}}
#' @param squash Plotting parameters. See \code{\link[wgaim]{link.map.cross}}
#' @param colqtl Color to plot QTL regions. See \code{\link[mpMap]{qtlmap}}
#' @param ... Additional arguments 
#' @return Modification of link.map.cross function from wgaim to allow more
#' general input objects and to highlight regions around QTL. If any markers
#' are labelled "QTLx" then they will be plotted in a different color.  
#' @seealso \code{\link[wgaim]{link.map.cross}}, \code{\link[mpMap]{qtlmap}}

plotlink.map <- 
function (object, chr, max.dist, marker.names = TRUE, tick = FALSE, 
    squash = TRUE, colqtl="red", ...) 
{
  circ <- function(x, y, shiftx = 0, shifty = 0, ely = 1, elx = 1)
    ((x - shiftx)^2)/elx + ((y - shifty)^2)/ely

    dots <- list(...)
    old.xpd <- par("xpd")
    par(xpd = TRUE)
    on.exit(par(xpd = old.xpd))
    if (inherits(object, "map")) map <- object else 
    if (inherits(object, "mpcross") & !is.null(object$map)) map <- object$map
    else stop("Incorrectly formatted object has been input\n")
#    map <- pull.map(object)
#    map <- object
    if (!missing(chr)) {
        if (any(is.na(pmatch(chr, names(map))))) 
            stop("Some names of chromosome(s) subset do not match names of map.")
        map <- map[chr]
    }
    n.chr <- length(map)
    mt <- list()
    if (!missing(max.dist)) 
        map <- lapply(map, function(el, max.dist) el[el < max.dist], 
            max.dist)
    maxlen <- max(unlist(lapply(map, max)))
    if (!marker.names) {
        chrpos <- 1:n.chr
        thelim <- range(chrpos) + c(-0.5, 0.5)
    }
    else {
      if (!is.na(pmatch("cex", names(dots)))) 
        cex <- dots$cex
      else cex <- par("cex")
      if(!squash)
        chrpos <- seq(1, n.chr * 3, by = 3)
      else
        chrpos <- seq(1, n.chr * 2, by = 2)
      thelim <- range(chrpos) + c(-1.6, 1.35)
      for (i in 1:n.chr) {
        mt[[i]] <- map[[i]]
        conv <- par("pin")[2]/maxlen
        for (j in 1:(length(mt[[i]]) - 1)) {
          ch <- mt[[i]][j + 1] * conv - (mt[[i]][j] * conv + 
                                         10 * par("csi") * cex/9)
          if (ch < 0) {
            temp <- mt[[i]][j + 1] * conv + abs(ch)
            mt[[i]][j + 1] <- temp/conv
          }
        }
      }
        maxlen <- max(unlist(lapply(mt, max)))
        names(mt) <- names(map)
    }
    plot(0, 0, type = "n", ylim = c(maxlen, 0), xlim = thelim, 
        xaxs = "i", ylab = "Location (cM)", xlab = "Chromosome", 
        axes = FALSE, ...)
    axis(side = 2, ylim = c(maxlen, 0))
    pins <- par()$plt
    for (i in 1:n.chr) {
        if (marker.names) {
            text(chrpos[i] + 0.5, mt[[i]], names(map[[i]]), adj = c(0, 
                0.5), ...)
	    m <- grep("QTL", names(map[[i]]))
            text(chrpos[i] + 0.5, mt[[i]][m], names(map[[i]])[m], adj = c(0, 
                0.5), col=colqtl, ...)
            segments(chrpos[i] + 0.25, map[[i]], chrpos[i] + 
                0.3, map[[i]])
            segments(chrpos[i] + 0.3, map[[i]], chrpos[i] + 0.4, 
                mt[[i]])
            segments(chrpos[i] + 0.4, mt[[i]], chrpos[i] + 0.45, 
                mt[[i]])
        }
        barl <- chrpos[i] - 0.03
        barr <- chrpos[i] + 0.03
        segments(barl, min(map[[i]]), barl, max(map[[i]]), lwd = 1)
        segments(barr, min(map[[i]]), barr, max(map[[i]]), lwd = 1)
        segments(barl - 0.17, map[[i]], barr + 0.17, map[[i]])
        xseq <- seq(barl, barr, length = 20) - chrpos[i]
        yseq <- circ(xseq, xseq, ely = 1, elx = 0.07/maxlen)
        yseq <- yseq - max(yseq)
        lines(xseq + chrpos[i], min(map[[i]]) + yseq)
        lines(xseq + chrpos[i], max(map[[i]]) - yseq)
    }
    axis(side = 1, at = chrpos, labels = names(map), tick = tick)
    if (is.na(pmatch("main", names(dots))) & !as.logical(sys.parent())) 
        title("Genetic Map")
    invisible(list(mt = mt, map = map, chrpos = chrpos))
}
