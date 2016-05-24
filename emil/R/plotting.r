#' Add vertical or horizontal lines to a plot
#'
#' @param x Coordinates of vertical lines.
#' @param lend Line ending style, see \code{\link{par}}.
#' @param ... Sent to \code{\link{segments}}.
#' @examples
#' plot(0:10, 0:10, type="n")
#' hlines(0:4*2.5, col="#dddddd")
#' points(0:10, 0:10)
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
vlines <- function(x, lend=1, ...)
    segments(x, par("usr")[3], x, par("usr")[4], lend=lend, ...)

#' @param y Coordinates of horizontal lines.
#' @rdname vlines
#' @export

hlines <- function(y, lend=1, ...)
    segments(par("usr")[1], y, par("usr")[2], y, lend=lend, ...)

#' Plots an axis the way an axis should be plotted.
#'
#' @param ... Sent to \code{\link{axis}}.
#' @param las Rotation of axis labels. Always horizontal by default.
#' @param lwd Width of the line drawn along the plot area. Omitted by default
#'   since it overlaps with \code{\link{box}} and causes it to look thicker
#'   where the axis is.
#' @param lwd.ticks Width of the tick lines. These are kept by default.
#' @param lend Line endings, see \code{\link{par}}.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
nice_axis <- function(..., las=1, lwd=0, lwd.ticks=par("lwd"), lend=2){
    axis(..., las=las, lwd=0, lwd.ticks=lwd.ticks, lend=1)
    if(lwd){
        # Plot the line along the axis
        args <- c(list(...), list(lwd=lwd, lwd.ticks=0, lend=lend))
        args$labels <- FALSE
        do.call(axis, args)
    }
}

#' Plots a box around a plot
#'
#' @param lend Line ending style, see \code{\link{par}}. Defaults to square.
#' @param ljoin Line joint style, see \code{\link{par}}. Defaults to mitre,
#'   i.e. 90 degree corners in this case.
#' @param ... Sent to \code{\link{box}}.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
nice_box <- function(lend=2, ljoin=1, ...) box(lend=lend, ljoin=ljoin, ...)

#' Get color palettes
#'
#' Can be used to modify an existing palette, e.g. change brightness,
#' or to generate a palette for a response vector.
#'
#' @param x Character vector of colors or factor of class memberships to
#'   generate colors for.
#' @return A character vector of hex colors.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
get_color <- function(x, ...){
    UseMethod("get_color")
}
#' @method get_color default
#' @param s Saturation. \code{s = 0} leaves it unchanged, \code{0 < s <= 1}
#'   increases, and \code{-1 <= s < 0} decreases.
#' @param v Value. \code{s = 0} leaves it unchanged, \code{0 < s <= 1}
#'   increases, and \code{-1 <= s < 0} decreases.
#' @param alpha Transparency.
#' @rdname get_color
#' @export
get_color.default <- function(x, s, v, alpha, ...){
    col <- as.data.frame(t(rgb2hsv(col2rgb(x))))
    if(!missing(alpha)) col$alpha <- alpha

    if(!missing(s)){
        if(length(s) != length(x) && length(s) != 1){
            warning("The length of `x` and `s` do not match, only using the first element.")
            s <- rep(s[1], length(x))
        }
        col$s <- ifelse(s < 0, col$s + s*col$s, col$s + s*(1-col$s))
    }
    if(!missing(v)){
        if(length(v) != length(x) && length(v) != 1){
            warning("The length of `x` and `v` do not match, only using the first element.")
            v <- rep(v[1], length(x))
        }
        col$v <- ifelse(v < 0, col$v + v*col$v, col$v + v*(1-col$v))
    }
    col[T] <- lapply(col, function(x) pmin(1, pmax(0, x)))
    
    structure(do.call(hsv, col), names=names(x))
}
#' @method get_color factor
#' @param levels If \code{TRUE} a palette with one color per level of \code{x}
#'   is returned. If \code{FALSE} one color per element in \code{x} is returned.
#' @param col Color palette with one color per class or the name of the color
#'   brewer palette to use, see \code{name} argument of \code{\link[RColorBrewer]{brewer.pal}}
#'   for a list of possible values.
#' @param ... Sent to \code{\link{get_color.default}}.
#' @rdname get_color
#' @export
get_color.factor <- function(x, levels=FALSE, col="Set1", ...){
    if(length(col) == 1){
        nice_require("RColorBrewer")
        suppressWarnings(col <- RColorBrewer::brewer.pal(1000, col))
    }
    if (length(levels(x)) > length(col)) 
        warning("Too few colors to assign unique ones to each class.")
    col <- rep(col, ceiling(length(levels(x))/length(col)))[seq_along(levels(x))]

    col <- get_color(col, ...)
    if(levels){
        structure(col, names=levels(x))
    } else {
        structure(col[as.integer(x)], names=names(x))
    }
}

