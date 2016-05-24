#' River Layout Coordinates Calculation
#' 
#' This function calculates best fit plotting coordinates to enable rivers to
#' be shown on river charts. The output is a list, which can be used when
#' plotting the river chart and the information on the river chart. It provides
#' an opportunity to change the coordinates and other plotting parameters
#' before actually plotting.
#' 
#' 
#' @param river a vector of river names.
#' @param length a vector of river lengths.
#' @param parent a vector of river parents. The parent of a river is the river
#' into which it flows.  The parent of the main stream is "NA".
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
#' (\code{direction = 1}), or from left to right (\code{direction = -1}). By
#' default, \code{direction = 1}.
#' @param margin a value. The margin height between rivers in the topological
#' plot. By default, \code{margin = 0.5} and margin height is 0.5 times the
#' river height in the river chart.
#' @return The \code{RiverLayout} returns a list containing the data for river
#' chart plotting. The list includes, \item{riverdata}{a data frame. This
#' contains input vectors \code{river}, \code{length}, \code{parent},
#' \code{position} and \code{distance}. It also includes calculated
#' x-coordinates of river mouths (\code{rmouth}) and sources (\code{rsource}),
#' defined in the same units as the inputs \code{length} and \code{distance}.
#' The last vector included is the row number for each river (\code{row}), in
#' which the main stream has a fixed row number of 0.} \item{H.MAX}{the number
#' of river rows.} \item{H.SIZE}{the height of each river row in the
#' topological plot.} \item{W.MAX}{the width of river layout, in the same units
#' as \code{length} and \code{distance}.} \item{W.SIZE}{the reciprocal of
#' \code{W.MAX}.} \item{X1}{normalised x-coordinate of river mouths.}
#' \item{X2}{normalised x-coordinate of river sources.} \item{Y}{normalised
#' y-coordinate of rivers.} \item{direction}{flow direction. Rivers flow from
#' right to left (\code{direction = 1}), or from left to right (\code{direction
#' = -1}).}
#' @note There is one and only one mainstream input for each function call.
#' @author Feng Mao
#' @seealso \code{\link{RiverDraw}}, \code{\link{RiverMap}}.
#' @keywords hplot
#' @examples
#' 
#' 
#' data(Ballinderry)
#' 
#' # River flows right
#' riverlayout <- RiverLayout(B.river$River,B.river$Length,
#'                            B.river$Parent,B.river$Position,
#'                            B.river$Distance, direction = -1)
#' 
#' # River flows left
#' riverlayout.left <- RiverLayout(B.river$River,B.river$Length,
#'                                 B.river$Parent,B.river$Position,
#'                                 B.river$Distance)
#' 
#' str(riverlayout)
#' 
#' @export RiverLayout
RiverLayout <- function(river, length, parent, position, distance,
                        row = NA,
                        direction = 1, 
                        margin = 0.5){
  
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
  
  list(riverdata = riverlayout, H.MAX = H.MAX, H.SIZE = H.SIZE, W.MAX = W.MAX, W.SIZE = W.SIZE, X1 = X1, X2 = X2, Y = Y, direction = direction)
  
}
