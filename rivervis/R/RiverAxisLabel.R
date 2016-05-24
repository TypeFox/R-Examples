#' River Axis Labels
#' 
#' This adds left or right axis labels to the river chart.
#' 
#' 
#' @param label the axis label to be shown on the river chart.
#' @param riverlayout the output list of \code{RiverLayout}.
#' @param cex text size.
#' @param adj text adjustment. One or two values in the range [0,1] for x and y
#' (optional) adjustment.
#' @param srt text angle.
#' @param col text colour.
#' @param pos text position. \code{1} for below, \code{2} for left, \code{3}
#' for above, and \code{4} for right. See \code{par} for details.
#' @param offset text position offset.
#' @param side left (\code{"L"}) or right (\code{"R"}) axis.
#' @param mainonly the axis title is only shown for the main stream only
#' (\code{"TRUE"}) or not (\code{"FALSE"}).
#' @author Feng Mao
#' @seealso \code{\link{RiverLayout}}, \code{\link{RiverDraw}},
#' \code{\link{RiverMap}}.  \code{\link{par}}.
#' @keywords hplot
#' @examples
#' 
#' 
#' # see examples below
#' data(Ballinderry)
#' 
#' riverlayout <- RiverLayout(B.river$River,B.river$Length,
#'                            B.river$Parent,B.river$Position,
#'                            B.river$Distance, direction = -1)
#' RiverDraw(riverlayout)
#' 
#' RiverBar(B.siteaspt$Site, B.siteaspt$River, B.siteaspt$Distance, 
#'          B.siteaspt[4:5], riverlayout, range = c(0,8), 
#'          bar.col = c("#5381FFFF", "#FF3931FF"), lbl.adj = c(0.5,1.3))
#' 
#' RiverTM(c(0,2,4,6,8,10), B.siteaspt[4:5], riverlayout, 
#'         pos=-1, side = "L", range = c(0,8), label = c(0,2,4,6,8))
#' 
#' RiverAxisLabel("ASPT score", riverlayout, adj = c(0.5, -3))
#' 
#' 
#' @export RiverAxisLabel
RiverAxisLabel <- function(label, riverlayout,
                           cex = 0.7,
                           adj = c(0.5, -2),
                           srt = 90,
                           col = "black",
                           pos = NULL,
                           offset = 0.5,
                           side = "L",
                           mainonly = TRUE){
  
  # Data transfer
  RIVER.DATA <- riverlayout[[1]]
  H.MAX <- riverlayout[[2]]
  H.SIZE <- riverlayout[[3]]
  W.MAX <- riverlayout[[4]]
  W.SIZE <- riverlayout[[5]]
  X1 <- riverlayout[[6]]
  X2 <- riverlayout[[7]]
  Y <- riverlayout[[8]]
  DIRECTION <- riverlayout[[9]]
  
  # left or right?
  
  if (DIRECTION == -1){
    X <- X1
    X1 <- X2
    X2 <- X 
  }
  
  
  # River Axis Title
  if (side == "L"){
    
    LOC.ROW <- data.frame(X.ROW = X1, Y.ROW = Y, ROW = RIVER.DATA$row)
    
    for (i in (1:length(Y))[duplicated(Y)]){
      LOC.ROW <- LOC.ROW[LOC.ROW$X.ROW != max(LOC.ROW$X.ROW[LOC.ROW$Y.ROW==Y[i]]) | LOC.ROW$Y.ROW!=Y[i],]
    }
    
    if (mainonly){
      text(LOC.ROW$X.ROW[LOC.ROW$ROW==0],
           LOC.ROW$Y.ROW[LOC.ROW$ROW==0] + 0.5 * H.SIZE,
           labels = label,
           cex = cex,
           adj = adj,
           col = col,
           srt = srt)
    }else{
      text(LOC.ROW$X.ROW,
           LOC.ROW$Y.ROW + 0.5 * H.SIZE,
           labels = label,
           cex = cex,
           adj = adj,
           col = col,
           srt = srt)
    }
    
  }
  
  if (side == "R"){
    
    LOC.ROW <- data.frame(X.ROW = X2, Y.ROW = Y, ROW = RIVER.DATA$row)
    
    
    for (i in (1:length(Y))[duplicated(Y)]){
      LOC.ROW <- LOC.ROW[LOC.ROW$X.ROW != max(LOC.ROW$X.ROW[LOC.ROW$Y.ROW==Y[i]]) | LOC.ROW$Y.ROW!=Y[i],]
    }
    
    if (mainonly){
      text(LOC.ROW$X.ROW[LOC.ROW$ROW==0],
           LOC.ROW$Y.ROW[LOC.ROW$ROW==0] + 0.5 * H.SIZE,
           labels = label,
           cex = cex,
           adj = adj,
           col = col,
           srt = srt)
    }else{
      text(LOC.ROW$X.ROW,
           LOC.ROW$Y.ROW + 0.5 * H.SIZE,
           labels = label,
           cex = cex,
           adj = adj,
           col = col,
           srt = srt)
    }
    
  }
}
