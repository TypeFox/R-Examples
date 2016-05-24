#' River Layout Coordinates Calculation and Plotting
#' 
#' This function calculates plotting coordinates for rivers and draws the river
#' chart accordingly.
#' 
#' 
#' @param river a vector of river names.
#' @param length a vector of river lengths.
#' @param parent a vector of river parents. The parent of a river is the river
#' into which it flows.  The parent of the main stream is NA.
#' @param position a vector of river positions. The river position indicates
#' its position relative to its parent - whether it is a left bank river, right
#' bank river or main stream. The left bank river is on the left when looking
#' downstream of its parent. The right bank river is on the right when looking
#' downstream of its parent.  The value of position can be "R", "L" or "M". For
#' the main stream, the value is "M".
#' @param distance a vector of distances denotes the distance between the
#' mouths of each river and the mouths of each river's parent.
#' @param row a vector of row numbers. The main stream is on row 0. In the
#' river chart, rivers with negative row numbers are plotted below the main
#' stream while rivers with positive row numbers are plotted above the main
#' stream. If a value for row is provided, the rivers will be plotted according
#' to the provided row numbers. If a value for row is not provided, a vector of
#' best fit row numbers will be calculated before plotting.
#' @param direction a value. In the river chart, rivers flow from right to left
#' (direction = 1), or from left to right (direction = -1). By default,
#' direction = 1.
#' @param margin a value. The margin height between rivers in the topological
#' plot. By default, margin = 0.5 and margin height is 0.5 times the river
#' height in the river chart.
#' @param bd.col river border colour.
#' @param bg.col background colour.
#' @param ln.col lead line colour.
#' @param ln.lty lead line style.
#' @param ln.lwd lead line width.
#' @param pt.shw show anchor point (\code{TRUE}) or not (\code{FALSE}). Anchor
#' points represent the locations of the river mouths.
#' @param pt.col anchor point colour.
#' @param pt.pch anchor point character.
#' @param pt.bg anchor point background(fill) colour when \code{pch=21:25}.
#' @param pt.cex anchor point size.
#' @param pt.lwd anchor point border width.
#' @param mar.t top margin size. This ranges in [0, 1] where 1 is the total
#' height of the diagram region.
#' @param mar.b bottom margin size.  This ranges in [0, 1] where 1 is the total
#' height of the diagram region.
#' @param mar.l left margin size.  This ranges in [0, 1] where 1 is the total
#' width of the diagram region.
#' @param mar.r right margin size. This ranges in [0, 1] where 1 is the total
#' width of the diagram region.
#' @return The \code{RiverMap} returns a list containing data for river map
#' drawing, and plots the river map accordingly. The output list can be used
#' for further plotting. The output list includes, \item{riverdata}{a data
#' frame. This contains input vectors \code{river}, \code{length},
#' \code{parent}, \code{position} and \code{distance}. It also includes
#' calculated x-coordinates of river mouths (\code{rmouth}) and sources
#' (\code{rsource}), defined in the same units as the inputs \code{length} and
#' \code{distance}. The last included vector is the row number for each river
#' (\code{row}), in which the main stream has a fixed row number of 0.}
#' \item{H.MAX}{the number of rows.} \item{H.SIZE}{the height of each row in
#' the topological plot.} \item{W.MAX}{the width of river layout, in the same
#' units as \code{length} and \code{distance}.} \item{W.SIZE}{the reciprocal of
#' \code{W.MAX}.} \item{X1}{normalised x-coordinate of river mouths.}
#' \item{X2}{normalised x-coordinate of river sources.} \item{Y}{normalised
#' y-coordinate of rivers.} \item{direction}{flow direction. Flow from right to
#' left (\code{direction = 1}), or from left to right (\code{direction = -1}).}
#' @author Feng Mao
#' @seealso \code{\link{RiverLayout}}, \code{\link{RiverDraw}},
#' \code{\link{par}}.
#' @keywords hplot
#' @examples
#' 
#' 
#' data(Ballinderry)
#' 
#' riverlayout <- RiverLayout(B.river$River,B.river$Length,
#'                            B.river$Parent,B.river$Position,
#'                            B.river$Distance, direction = -1)
#' str(riverlayout)
#' 
#' RiverMap(B.river$River,B.river$Length,B.river$Parent,
#'          B.river$Position, B.river$Distance)[[1]]
#' 
#' RiverMap(B.river$River,B.river$Length,B.river$Parent,
#'          B.river$Position, B.river$Distance, 
#'          row = c(5,-1,6,3,-4,2,-6,7), direction = -1)
#' 
#' @export RiverMap
RiverMap <- function(river, length, parent, position, distance, 
                     row = NA, 
                     direction = 1, # when Direction=1, river mouth on left. when Direction=-1, river mouth on right.
                     margin = 0.5,
                     bd.col = "black",
                     ln.col = "grey40",
                     ln.lty = 3,
                     ln.lwd = 1,
                     bg.col = "grey80",
                     pt.shw = TRUE,
                     pt.col = "black",
                     pt.pch = 20, 
                     pt.bg = "black",
                     pt.cex = 1,
                     pt.lwd = 1,
                     mar.t = 0.05,  
                     mar.b = 0.05,
                     mar.l = 0.2,
                     mar.r = 0.1
){ 
  
  par(mar=c(0,0,0,0))
  
  riverlayout <- data.frame(river = river, length = length, parent = parent, position = position, distance = distance, stringsAsFactors = FALSE)
  
  OBN <- nrow(riverlayout) # Observation number, or number or rows
  
  riverlayout <- cbind(rivercode = paste("river",c(1:OBN),sep=""), riverlayout, stringsAsFactors = FALSE) # Allocate rivercode for rivers
  
  riverlayout$parent[riverlayout$position == "M"] <- NA # Make sure the Parent of mainstream is NA
  
  path <- PathBuild(riverlayout$river, riverlayout$parent, OBN)
  
  DIGITMAX <- ncol(path)-1
  
  pos <- RelPos(path, riverlayout, OBN, DIGITMAX)
  
  digitweight <- DigitWeight(DIGITMAX)
  
  posmatrix <- RelPosMatrix(pos, DIGITMAX)
  
  # Calculate Mouth and Source
  
  riverlayout <- cbind(riverlayout, MouthSource(path, riverlayout, OBN))
  
  # Calculate Row
  
  
  if(all(is.na(row))){
    riverlayout <- merge(riverlayout, RowCal(posmatrix, digitweight, riverlayout, path, OBN, DIGITMAX), by = "river", sort = FALSE)
    
    row <- riverlayout$row
    
  } else{  
    riverlayout <- data.frame(riverlayout, row = row)
  }
  
  row <- riverlayout$row
  rsource <- riverlayout$rsource
  rmouth <- riverlayout$rmouth
  
  # Judge flow direction
  if (direction == -1){
    MAX.SOURCE <- max(rsource)
    row <- -row
    rsource <- MAX.SOURCE - rsource
    rmouth <- MAX.SOURCE - rmouth
    riverlayout$rsource <- rsource
    riverlayout$rmouth <- rmouth
    riverlayout$row <- row
  }
  
  # Calculate unit height
  H.MAX <- max(row)-min(row)+1 # total number of rows
  
  H.SIZE <- 1/(H.MAX + H.MAX*margin + 1) # define the unit height of each river. Assume the margin between rows is HSIZE/2
  
  Y.ZERO <- abs(min(row[row<=0])) * (margin*H.SIZE + H.SIZE) + margin*H.SIZE # get the y coordinate for row 0 as a reference line
  
  Y <- row * (margin*H.SIZE + H.SIZE) + Y.ZERO # Y of left bottom points of river rectangles
  
  # Calculate unit width
  W.MAX <- max(rsource, rmouth) - min(rsource, rmouth) # maximum width in original units (km)
  
  W.SIZE <- 1/W.MAX # leave some space for the right
  
  X1 <- rmouth * W.SIZE # X of Mouth location
  
  X2 <- rsource * W.SIZE # X of Source location
  
  # Plot new sheet
  plot.new()
  
  # plotting margin
  par(usr = c(-mar.l, 1+mar.r, -mar.b, 1+mar.t))
  
  # Plot river rectangles
  rect(X1, Y, X2, Y+H.SIZE, col = bg.col, border = bd.col) # draw river rectangles
  
  # Plot lead line
  Y.PARENT <- Y[match(riverlayout$parent, riverlayout$river)] # Y of Parent of each river, using dictionary technique
  
  segments(X1, Y, X1, Y.PARENT, col = ln.col, lty = ln.lty, lwd = ln.lwd)
  
  # Plot river rectangle-frames
  rect(X1, Y, X2, Y+H.SIZE, col = NA, border = bd.col) # draw the frame of river rectangles, in case they have been covered by leadlines
  
  # Plot anchor points
  if (pt.shw){
    points(X1, Y.PARENT, pch=pt.pch, bg = pt.bg, cex = pt.cex, lwd = pt.lwd, col = pt.col) # plot anchor points
  }
  
  list(riverdata = riverlayout, H.MAX = H.MAX, H.SIZE = H.SIZE, W.MAX = W.MAX, W.SIZE = W.SIZE, X1 = X1, X2 = X2, Y = Y, direction = direction)
  
}
