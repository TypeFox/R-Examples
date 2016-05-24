#' Tick Marks on River Charts
#' 
#' This adds tick marks to the river chart.
#' 
#' 
#' @param tickmark a vector of tick mark values.
#' @param value the variables which the tick marks are for.
#' @param riverlayout the output list of \code{RiverLayout} or \code{RiverMap}.
#' @param range bar-chart value range. A vector of two values indicating lower
#' limit and upper limit.
#' @param lbl.shw show labels of tick marks (\code{TRUE}) or not
#' (\code{FALSE}).
#' @param lbl.col label colour.
#' @param lbl.cex label size.
#' @param lbl.row show one label per row (\code{TRUE}) or not (\code{FALSE}).
#' @param label a vector of tick mark labels.
#' @param side position of tick marks. \code{"l"} for left and \code{"r"} for
#' right.
#' @param pos position of tick marks. \code{-1} for in and \code{1} for out.
#' @param tm.l tick mark length.
#' @param tm.col tick mark colour.
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
#' RiverPoint(NA,B.elevation$River, B.elevation$Distance, 
#'            B.elevation$Elevation, riverlayout)
#' 
#' RiverTM(c(0, 100, 200, 300, 400, 500), B.elevation[3], riverlayout, 
#'         pos=-1, side = "R", range = c(0,500), 
#'         label = c(0, 100, 200, 300, 400, 500))
#' 
#' @export RiverTM
RiverTM <- function(tickmark, # a vector of tick mark values
                    value, # original data
                    riverlayout,
                    range = NA,
                    side = "L", # tick mark on left or right
                    pos = 1, # in (-1) or out (1)
                    tm.l = 1,
                    tm.col = "black",
                    lbl.shw = TRUE,
                    lbl.col = "black",
                    lbl.cex = 0.7,
                    lbl.row = TRUE,
                    label = NA){ # relative length of the tick mark
  
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
  
  #TickMark
  if (all(is.na(range))){
    tickmark <- tickmark[which(tickmark>(min(value)-(max(value)-min(value))*0.05/0.9) & tickmark<(max(value)+(max(value)-min(value))*0.05/0.9))]
  } else{
    tickmark <- tickmark[which(tickmark>(min(range)-(max(range)-min(range))*0.05/0.9) & tickmark<(max(range)+(max(range)-min(range))*0.05/0.9))]
  }
  
  
  L.TM <- 0.005 * tm.l # the default length of tick mark is 0.005
  
  # left or right?
  
  if (DIRECTION == -1){
    X <- X1
    X1 <- X2
    X2 <- X 
  }
  
  
  if(all(is.na(label))){
    
    label <- tickmark
  }
  
  
  if (all(is.na(range))){
    VALUE.MAX <- max(value)
    VALUE.MIN <- min(value)
  } else{
    VALUE.MAX <- max(range)
    VALUE.MIN <- min(range)
  }
  
  VALUE.SIZE <- H.SIZE * 0.9/(VALUE.MAX-VALUE.MIN)
  
  if (side == "L"){
    
    LOC.ROW <- data.frame(X.ROW = X1, Y.ROW = Y)
    
    if (lbl.row){
      for (i in (1:length(Y))[duplicated(Y)]){
        LOC.ROW <- LOC.ROW[LOC.ROW$X.ROW != max(LOC.ROW$X.ROW[LOC.ROW$Y.ROW==Y[i]]) | LOC.ROW$Y.ROW!=Y[i],]
      }
    }
    
    X.TM.LBL <- rep(LOC.ROW$X.ROW, each=length(tickmark))
    
    Y.TM.LBL <- rep(LOC.ROW$Y.ROW, each=length(tickmark)) + rep((tickmark - VALUE.MIN) * VALUE.SIZE, nrow(LOC.ROW)) + 0.05 * H.SIZE
    
    X.TM <- rep(X1, each=length(tickmark))
    
    Y.TM <- rep(Y, each=length(tickmark)) + rep((tickmark - VALUE.MIN) * VALUE.SIZE, nrow(RIVER.DATA)) + 0.05 * H.SIZE
    
    segments(X.TM, 
             Y.TM, 
             X.TM - pos * L.TM, 
             Y.TM, col = tm.col)
    
    if (lbl.shw){
      text(X.TM.LBL - L.TM,
           Y.TM.LBL,
           labels = label,
           cex = lbl.cex,
           adj = c(1,0.5),
           col = lbl.col)
    }
    
  }
  
  if (side == "R"){
    
    LOC.ROW <- data.frame(X.ROW = X2, Y.ROW = Y)
    
    if (lbl.row){
      for (i in (1:length(Y))[duplicated(Y)]){
        LOC.ROW <- LOC.ROW[LOC.ROW$X.ROW != min(LOC.ROW$X.ROW[LOC.ROW$Y.ROW==Y[i]])| LOC.ROW$Y.ROW!=Y[i],]
      }
    }
    
    X.TM.LBL <- rep(LOC.ROW$X.ROW, each=length(tickmark))
    
    Y.TM.LBL <- rep(LOC.ROW$Y.ROW, each=length(tickmark)) + rep((tickmark - VALUE.MIN) * VALUE.SIZE, nrow(LOC.ROW)) + 0.05 * H.SIZE
    
    X.TM <- rep(X2, each=length(tickmark))
    
    Y.TM <- rep(Y, each=length(tickmark)) + rep((tickmark - VALUE.MIN) * VALUE.SIZE, nrow(RIVER.DATA)) + 0.05 * H.SIZE
    
    segments(X.TM, 
             Y.TM, 
             X.TM + pos * L.TM, 
             Y.TM, col = tm.col)
    
    if (lbl.shw){
      text(X.TM.LBL + L.TM,
           Y.TM.LBL,
           labels = label,
           cex = lbl.cex,
           adj = c(0,0.5),
           col = lbl.col)
    }
    
  }
  
  
}
