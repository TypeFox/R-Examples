#' River Labels on River Charts
#' 
#' This adds the name labels to the plotted rivers.
#' 
#' 
#' @param riverlayout the output list of \code{RiverLayout} or \code{RiverMap}.
#' @param corner river label position, which can be at any of the four river
#' chart corners. \code{"lt"} for left-top, \code{"lb"} for left-bottom,
#' \code{"rt"} for right-top, \code{"rb"} for right-bottom.
#' @param cex text size.
#' @param adj text adjustment. One or two values in [0,1] for x and y
#' (optional) adjustment.
#' @param srt text angle.
#' @param col text colour.
#' @param pos text position. \code{1} for below, \code{2} for left, \code{3}
#' for above, and \code{4} for right. See \code{par} for details. It overrides
#' \code{adj} if it is not \code{NULL}.
#' @param offset text position offset.
#' @author Feng Mao
#' @seealso \code{\link{RiverLayout}}, \code{\link{RiverDraw}},
#' \code{\link{RiverMap}}.  \code{\link{par}}.
#' @keywords hplot
#' @examples
#' 
#' 
#' data(Ballinderry)
#' 
#' riverlayout <- RiverLayout(B.river$River,B.river$Length,
#'                            B.river$Parent,B.river$Position,
#'                            B.river$Distance, direction = -1)
#' RiverDraw(riverlayout)
#' 
#' RiverLabel(riverlayout, corner = "lt", srt = 0, adj = c(0, -0.7))
#' 
#' RiverLabel(riverlayout, corner = "lb")
#' 
#' RiverLabel(riverlayout, corner = "rt", srt = -90)
#' 
#' @export RiverLabel
RiverLabel <- function(riverlayout,
                       cex = 0.7,
                       adj = c(0, -1),
                       srt = 90,
                       col = "black",
                       pos = NULL,
                       offset = 0.5,
                       corner = "lb"     # left-top = lt, left-bottom = lb, right-top = rt, right-bottom = rb
){ 
  
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
  
  if (DIRECTION == -1){
    if (corner == "lb"){
      X.RIVER <- X2
      Y.RIVER <- Y
    }
    if (corner == "lt"){
      X.RIVER <- X2
      Y.RIVER <- Y + H.SIZE
    }
    if (corner == "rb"){
      X.RIVER <- X1
      Y.RIVER <- Y
    }
    if (corner == "rt"){
      X.RIVER <- X1
      Y.RIVER <- Y + H.SIZE
    }
    
  }else{
    if (corner == "lb"){
      X.RIVER <- X1
      Y.RIVER <- Y
    }
    if (corner == "lt"){
      X.RIVER <- X1
      Y.RIVER <- Y + H.SIZE
    }
    if (corner == "rb"){
      X.RIVER <- X2
      Y.RIVER <- Y
    }
    if (corner == "rt"){
      X.RIVER <- X2
      Y.RIVER <- Y + H.SIZE
    }
    
  }
  
  text(X.RIVER, Y.RIVER, labels = RIVER.DATA$river, cex = cex, adj = adj, srt = srt, col = col, pos = pos, offset = offset)
  
}
