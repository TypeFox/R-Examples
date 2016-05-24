#' plot multcomp graphics
#' 
#' Plot graphic(s) for multcompTs or multcompLetters objects
#' 
#' The requested graphic is either plotted by itself or added to an existing
#' plot as specified by the arguments.  The placement can be controlled by
#' 'fig' and 'at'.
#' 
#' The fontsize and fontface of a plot of a multcompLetters object with 'type =
#' "Letters"' can be adjusted as describe on the \code{\link[grid]{gpar}} help
#' page.
#' 
#' @param x an object of class 'multcompTs' or 'multcompLetters'.
#' @param horizontal A logical scalar indicating whether the list of items
#' compared reads left to right (horizontal = TRUE) or top to bottom
#' (horizontal = FALSE).  If this multcomp graphic accompanies boxplots for
#' different levels or groups compared, the 'boxplot' argument 'horizontal' is
#' the negation of the multcomp plot 'horizontal' argument.
#' @param col The color for each group of items or factor levels.  The colors
#' will cross the different items or factor levels and will therefore have the
#' orientation specified via 'horizontal'.  If the number of columns exceeds
#' length(col), col is recycled.  For alternative choices for col, see "Color
#' Specification" in the \code{\link{par}} help page.
#' @param type An alternative display for either multcompTs or multcompLetters
#' is 'boxes' (or rectangles).  If type="boxes" with "multcompTs", the
#' "base(s)" of each "T" will be indicated by a triangle.
#' @param orientation The "standard" orientation has the 'multcompTs' pointing
#' towards the names of the items or factor levels; with the "reverse"
#' orientation, the bases of the "Ts" point away.  By default, the names are on
#' the left or below unless the mean of the relevant 'fig' range is less than
#' 0.5.
#' @param add TRUE to add to an existing plot; FALSE to start a new plot.  The
#' names of the factor levels or items compared will be plotted only if
#' add=FALSE.
#' @param at A numeric vector or matrix or a list with components "at" and
#' "width".  If a list, both components must be either a numeric vector or
#' matrix.  The numeric vector "at" (whether the function argument or "at"
#' component of the "at" list) must be either a numeric vector or matrix giving
#' the locations where the "Ts" or "Letters" graphics should be drawn.
#' length(at) is 1, 2 or 3 times the number of the number of factor levels or
#' items compared.  If length(at) is twice the number of levels or items
#' compared, it gives the range of the display for that level; the base of a
#' "T" will be at the midpoint.  If length(at) is three times the number of
#' items compared, the intermediate number will be the center of the base of
#' the "T".
#' @param width A numeric vector or matrix with as many rows as "Ts" or
#' "groups" and with up to three columns.  With one column, it will be the
#' "center" of the plot range for that group.  With two columns, they will
#' delimit the range.  With three, they will provide "bottom", "center", and
#' "top" of the range for that set of grouping indicators.  If "at" is a list,
#' the argument "width" is ignored and is taken from the list "at".
#' @param fig A numerical vector of the form 'c(x1, x2, y1, y2)' giving the
#' coordinates of the figure region in the display region of the plot device,
#' as described on the \code{\link{par}} help page.
#' @param lwd width of line to connect elements of "T" graphics that might not
#' otherwise be connected.
#' @param label.levels NA for no labels or distance away from the plot for the
#' labels as a proportion of the plot range.
#' @param label.groups NA for no labels or distance away from the plot for the
#' labels as a proportion of the plot range.
#' @param T.base A numeric scalar giving the proportion of the available space
#' devoted to the base of the Ts; used only when type="Ts".
#' @param ...  graphical parameters can be given as described on the
#' \code{\link{plot}} help page or for plot.multcompLetters as describe on the
#' \code{\link[grid]{gpar}} help page.
#' @return A list with two components: \item{at }{ A matrix with three columns
#' giving the middle and extremes of the display for each of the factor levels
#' or items compared.  } \item{width }{ A matrix with as many rows as "Ts" or
#' comparitor levels and with two columns giving the plot range for that
#' comparitor level.  }
#' @author Spencer Graves
#' @seealso \code{\link{multcompTs}} \code{\link{multcompLetters}}
#' \code{\link{multcompBoxplot}} \code{\link[grid]{gpar}}
#' @keywords aplot
#' @examples
#' 
#' ##
#' ## plot(multcompTs(...))
#' ##
#' dif4 <- c(.1, .02, .03, 1)
#' names(dif4) <- c("A-B", "A-C", "B-C", "A-D")
#' (mcT4 <- multcompTs(dif4))
#' # Standard plot, base of "Ts" point left
#' mcT4.1 <- plot(mcT4, label.groups=0.05)
#' # Redo using "at" = list 
#' plot(mcT4, label.groups=0.05, at=mcT4.1)
#' # Same plot with group labels closer to the figure
#' plot(mcT4, label.groups=0.02)
#' 
#' # Base of "Ts" point right 
#' plot(mcT4, label.groups=TRUE, orientation="r")
#' # Base of "Ts" point down 
#' plot(mcT4, horizontal=TRUE, label.groups=0.05)
#' # Base of "Ts" point up
#' plot(mcT4, horizontal=TRUE, label.groups=0.05,
#'      orientation="r")
#' 
#' # Same 4 plots but with boxes & triangles, not Ts
#' plot(mcT4, label.groups=0.05, type="b")
#' plot(mcT4, label.groups=0.05, orientation="r",
#'      type="b")
#' plot(mcT4, horizontal=TRUE, label.groups=0.05,
#'      type="b")
#' plot(mcT4, horizontal=TRUE, label.groups=0.05,
#'      orientation="r", type="b")
#' 
#' ##
#' ## plot(multcompLetters(...))
#' ##
#' # ... using dif4 from above
#' (mcL4 <- multcompLetters(dif4, Letters=LETTERS))
#' # Standard plot
#' \dontrun{
#' # Requires (grid)
#' mcL4.1 <- plot(mcL4, label.groups=0.05)
#' # Redo using "at" = list 
#' plot(mcL4, label.groups=0.05, at=mcL4.1)
#' 
#' # With bold face and larger font 
#' plot(mcL4, label.groups=0.05,
#'      fontsize=28, fontface="bold")
#' 
#' # Horizontal rather than vertical 
#' plot(mcL4, horizontal=TRUE, label.groups=0.05)
#' }
#' 
#' # Same as boxes rather than letters 
#' plot(mcL4, label.groups=0.05, type="b")
#' plot(mcL4, horizontal=TRUE, label.groups=0.05,
#'      type="b")
#' @importFrom graphics par text lines polygon rect
#' @export

"plot.multcompTs" <-
function(x,
   horizontal=FALSE, col=1:6,
   type=c("Ts", "boxes"), 
   orientation=c("standard", "reverse"),
   add=FALSE, at, width, fig=c(0, 1, 0, 1),
   lwd=3, label.levels=if(add)NA else 0.05,
   label.groups=NA, T.base=0.4, ...){
##
## 1.  Get row and column names
##
  obj <- deparse(substitute(x))
  lvls <- dimnames(x)[[1]]
  gps <- dimnames(x)[[2]]
  n <- length(lvls)
  k <- length(gps)
  if((n==0) | (k==0))
    stop("dimnames required for ", obj)
  col <- rep(col, len=k)
##
## 2.  Get plotting positions for lvls and gps
##
  At <- function(x, lvls, d){
    N <- length(lvls)
    if(length(x)<=N)      
      x <- t(outer(c(-0.4, 0, 0.4), x, "+"))
    else{
      x <- apply(x, 2, sort)
      if(length(x)<(3*N))
        x <- cbind(x[, 1], mean(x, 2, mean),
                   x[, 2])
    }
    dimnames(x) <- list(lvls,
         c("bottom", "center", "top"))
    x
  }
#
  at.list <- FALSE
  missW <- missing(width)
  {
    if(missing(at)){
      at <-{
        if(horizontal) 1:n else n:1
      }
    }
    else 
      if(is.list(at)){
        at.list <- TRUE
        if(is.null(names(at)))
          names(at) <- c("at", "width")
        width <- at$width
        at <- at$at
      }
  }
#
  if(!at.list){
    if(missW){
      width <- {
        if(horizontal) k:1 else 1:k
      }
    }      
    width <- At(width, gps)
    at <- At(at, lvls)
  }
  or. <- match.arg(orientation)
## 
## 3.  Set up the plot 
##
  op <- {
    if(add)par(fig=fig, xpd=NA, new=TRUE)
    else par(fig=fig, xpd=NA)
  }
  on.exit(par(op))
#
  {
    if(match.arg(type)=="Ts")
      plotTs(obj=x, at=at, width=width, 
         horizontal=horizontal, col=col,
         add=add, lwd=lwd,
         label.levels=label.levels, 
         label.groups=label.groups,
         T.base=T.base, orientation=or.,
         ...)
    else
      plotBoxes(obj=x, at=at, width=width, 
         horizontal=horizontal, col=col,
         add=add, label.levels=label.levels, 
         label.groups=label.groups,
         orientation=or., ...)
  }
  list(at=at, width=width)
}

