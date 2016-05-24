#' plot multcomp displays
#' 
#' Helper functions for plot.multcompTs and plot.multcompLetters.  These not
#' intended to be called directly and are hidden in a namespace.  You can use
#' 'getAnywhere' to see them.
#' 
#' The requested graphic is either plotted by itself or added to an existing
#' plot as specified by the arguments.
#' 
#' "plotTs" and "plotBoxes" use traditional R graphics and will not be
#' discussed further here.
#' 
#' "plotLetters" uses 'grid' graphics, because it seems to provide more support
#' for controlling the side-by-side placement of "Letters" of possibly
#' different colors and widths.  The "Letters" display will be positioned in
#' the "plot region" defined by fig and mar, assuming the entire device region
#' is 37 lines both wide and tall.  Thus, the plot region is diff(fig[1:2])*37
#' lines wide and diff(fig(1:2])*37 lines high.  If, for example, fig = c(0.9,
#' 1, 0, 1), this makes the plot region 3.7 lines wide.  With the default
#' mar=c(5, 4, 4, 2)+0.1 lines, the "width" of the plot region is therefore 3.7
#' - (4.1+2.1) = (-2.5) lines.  "plotLetters" initially ignores this
#' contradictory negative width, and centers the plot at the midpoint of h0 =
#' fig[1]+mar[2]/37, h1 = fig[2]-mar[4]/37, v0 = fig[3]+mar[1]/37, and v1 =
#' fig[4]-mar[3]/37.  The "Letters" for the different levels compared are
#' rescaled from at[, "center"] to fit inside At.rng = if(horizontal) c(h0, h1)
#' else c(v0, v1).  With "n" levels compared and at.rng = range(at[,
#' "center"]), at[, "center"] is expanded to (at.rng+/-0.5) and rescaled to
#' match At.rng; if(diff(At.rng)<=0), an error message is issued.
#' 
#' Meanwhile, the "Letters" are centered at the midpoint of W.rng =
#' if(horizontal) c(v0, v1) else v(h0, h1) [the opposite of At.rng]; the
#' argument "width" used by plotTs and plotBoxes is not used (and not even
#' accepted) by plotLetters.  If(label.levels), these are positioned in the
#' midpoint of the right margin in the "W" direction.
#' 
#' @aliases plotLetters
#' @param obj a matrix describing which levels (rows) will be plotted with
#' which groups (columns).  For plotTs and plotBoxes, obj is a matrix of
#' numbers from (-1, 0, 1).  For plotLetters, obj is a logical matrix = TRUE if
#' that "letter" (group or column of obj) is to be plotted with that level (row
#' of obj).
#' @param at an array with one row for each level and 3 columns giving low,
#' middle and high levels for the display for that level.
#' @param width an array with one row for each group of levels in the display
#' and 3 columns giving low, middle and high levels for the display for that
#' group.
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
#' @param add TRUE to add to an existing plot; FALSE to start a new plot.  The
#' names of the factor levels or items compared will be plotted only if
#' add=FALSE.
#' @param lwd line width for the display outline.
#' @param label.levels Distance from the plot region to print the names of the
#' levels as a proportion of the plot range; NA for no level labels.
#' @param label.groups Distance from the plot region to print the names of the
#' groups as a proportion of the plot range; NA for no level labels.
#' @param T.base A numeric scalar giving the proportion of the available space
#' devoted to the base of the Ts.
#' @param orientation If 'reversed', the base(s) of each "T" or traingle
#' indicating the master level(s) of that "undifferentiated class" will point
#' right or up (depending on horizontal) rather than down or left.
#' @param font.family character string naming the font family used by
#' "plotLetters".  This function plots the different "Letters" in different
#' colors by plotting one color at a time.  It's currently not smart enough to
#' align the letters properly except by assuming a mono-spaced font.
#' @param fig figure region = (x0, x1, y0, y1) as a proportion of the device
#' region.
#' @param mar margin = (lower, left, upper, right) in lines.
#' @param ...  graphical parameters can be given as described on the
#' \code{\link{plot}} help page or the \code{\link[grid]{gpar}} help page.
#' @return "Done"
#' @author Spencer Graves
#' @seealso \code{\link{plot.multcompTs}} \code{\link{plot.multcompLetters}}
#' \code{\link[grid]{gpar}}
#' @keywords aplot internal
#' @examples
#' 
#' # Designed to be called from plot.multcompTs
#' # or plot.multcompLetters, NOT directly by users.  
#' @importFrom grid grid.newpage viewport pushViewport unit textGrob gpar grobHeight
#' @importFrom grid convertHeight grobWidth convertWidth grid.text popViewport
#' @export

"plotLetters" <- 
function(obj, at, horizontal, col, add,
    label.levels, font.family="mono",
    fig=par("fig"), mar=par("mar"), ...){
##
## 1.  Modify "par"
##
  if(!add) grid::grid.newpage()
#  
  op <- par(family=font.family)
  on.exit(par(op))  
##
## 2.  Set up 'grid' graphics 
##  
# Use 'grid' graphics for 2 reasons:
#   (1) I couldn't figure out how to get 
#       consistent alignment with overlays
#       using 'text'
#   (2) 'grid' will compute the size 
#       of characters, so I can use that
#       to determine plotting positions.  
  vpFig <- grid::viewport(name="vp.fig")
grid::pushViewport(vpFig)
#  
  angle <- 90*horizontal
  lvl.rng <- range(at[, "center"])
##
## 3.  Compute figure and plot regions
##
  n <- dim(obj)[1]
  k <- dim(obj)[2]
# Figure region 
  W.fig <- fig[(1:2)+2*horizontal]
  At.fig <- fig[(3:4)-2*horizontal]
# Figure margins
  At.mar <- mar[c(1, 3)+horizontal]
  W.mar <- mar[c(2, 4)-horizontal]
# Plot region
  At.rng <- (At.fig+c(1, -1)*At.mar/37)
  W.rng <- (W.fig+c(1, -1)*W.mar/37)
##
## 4.  Rescale "at" to Grid's "npc" 
##     = "normalized parent coordinates"
##     within fig and mar
##  
  at.rng <- (range(at[, "center"])+c(-0.5, 0.5))
  d.at <- diff(at.rng)
  d.At <- diff(At.rng)
  At <- (At.rng[1]+d.At*(
              at[, "center"]-at.rng[1])/d.at)
#  At.npc <- unit(At, "npc")
#  At.nat <- unit(At, "native")
##
## 5.  Convert "Letters" to grobs =
##     Grid "graphics objects" 
##  
  n <- dim(obj)[1]
  k <- dim(obj)[2]  
  Ltrs <- dimnames(obj)[[2]]
  k1 <- k+1
  LtrsM <- c(Ltrs, "M")
  LtrsM. <- c(Ltrs, "Ref.M.")
  gLtrs <- vector(k1, mode="list")
  names(gLtrs) <- LtrsM
#  npc0 <- unit(0, "npc")
  nat0 <- grid::unit(0, "native")
  for(j in 1:k1)
    gLtrs[[j]] <- grid::textGrob(LtrsM[j], rot=angle,x=nat0,
       y=nat0, name=LtrsM.[j], gp=grid::gpar(col=col[j], ...))
#          name=LtrsM.[j], gp=gpar(...))
##
## 6.  Compute character widths
##
  wLtrs <- rep(NA, k1)
  names(wLtrs) <- LtrsM
  {
    if(horizontal)
      for(j in 1:k1){
        gH.j <- grid::grobHeight(gLtrs[[j]])
        wLtrs[j] <- grid::convertHeight(gH.j, "native",
                       valueOnly=TRUE)
      }
    else
      for(j in 1:k1){
        gW.j <- grid::grobWidth(gLtrs[[j]])
        wLtrs[j] <- grid::convertWidth(gW.j, "npc",
                       valueOnly=TRUE)
      }
  }
##
## 7.  Rescale to W.rng
##
#     7.1.  maxX = max width including Ref.M.       
  maxW <- max(wLtrs)
#     7.2.  w.Ltrs = adj. width excl. Ref.M.      
  w.Ltrs <- 0.5*(wLtrs[-k1]+maxW)
#     7.3.  Convert to a scale in "npc"      
  sumW <- cumsum(w.Ltrs)
#     7.4.  width(in "npc") = adjustment to W.rng
  w0 <- mean(W.rng)
  W <- (w0+sumW-mean(range(sumW)))
#  W.npc <- unit(W, "npc")
#  W.nat <- unit(W, "native")
##
## 8.  Plot
##  
  {
    if(horizontal)
      for(j in 1:k){
        n.j <- sum(obj[, j])
        Ltr.j <- rep(Ltrs[j], n.j)
        W.j <- rep(W[j], n.j)
        At.j <- At[obj[, j]]
        grid::grid.text(Ltr.j, At.j, W.j, rot=90,
                  gp=grid::gpar(col=col[j], ...))
      }
    else
      for(j in 1:k){
        n.j <- sum(obj[, j])
        Ltr.j <- rep(Ltrs[j], n.j)
        W.j <- rep(W[j], n.j)
        At.j <- At[obj[, j]]
        grid::grid.text(Ltr.j, W.j, At.j,
                  gp=grid::gpar(col=col[j], ...))
      }
  }
##
## 9.  Label the levels?    
##
  if(!is.na(label.levels)){
    lvls <- dimnames(obj)[[1]]
#   W.mar = (W.rng[1]- 0.5*dW - label.levels)    
    W.mar <- (W[1]-0.5*wLtrs[1]-label.levels)
    W.n <- rep(W.mar, n)
    if(horizontal){
      grid::grid.text(lvls, At, W.n, rot=90)
    }else{
      grid::grid.text(lvls, W.n, At)
    }
  }
##
## 10.  Clean up and quit.
##      
grid::popViewport()
#  
  "Done"
}

