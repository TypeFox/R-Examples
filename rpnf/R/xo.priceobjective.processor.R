#' This function adds Vertical Price Objectives calculated with the
#' Bullish Breakout and Bearish Breakdown Method (BM) to an P&F Table.
#' 
#' Finding the appropriate price objectives has been explained very good at 
#' http://stockcharts.com/school/doku.php?id=chart_school:chart_analysis:point_and_figure_pri, 
#' but this documentation is no longer available.
#' The function adds columns vpo_bm_boxnumber and vpo_bm_price to the given
#' P&F Table. vpo_bm_bonumber contains the boxnumber of the price objective,
#' while vpo_bm_price contains the real price objective.
#' 
#' @param data Input data
#' @param reversal Number of boxes for reversal
#' @param boxsize Size of one box
#' @param log Use logarithmic scale
#'
xo.priceobjective.processor <- function(data,reversal,boxsize,log) {
  # add new column to store vertical price objective for breakout method (BM)
  data$vpo_bm_boxnumber <- rep(NA,times=nrow(data))
  data$vpo_bm_price <- rep(NA,times=nrow(data))
  # add new column to store vertical price objective for reversal method (RM)
  data$vpo_rm_boxnumber <- rep(NA,times=nrow(data))
  data$vpo_rm_price <- rep(NA,times=nrow(data))
  # loop over every in set
  for (i in 1:nrow(data)) {
    # subset current data
    # FIXME this is a bottleneck
    mydata <- data[1:i,]
    # determine vertical price objective with breakout method
    # FIXME this is a bottleneck
    vert.obj <-  currentVPOBreakoutMethod(data=mydata,
                                           reversal=reversal,
                                           boxsize=boxsize,
                                           log=log) 
    data$vpo_bm_boxnumber[i] <- vert.obj$boxnumber
    data$vpo_bm_price[i] <- vert.obj$price
    # determine vertical price objective with breakout method
    # FIXME this is a bottleneck
    vert.obj <-  currentVPOReversalMethod(data=mydata,
                                           reversal=reversal,
                                           boxsize=boxsize,
                                           log=log) 
    data$vpo_rm_boxnumber[i] <- vert.obj$boxnumber
    data$vpo_rm_price[i] <- vert.obj$price
  }
  
  data
  # TODO implement function
}
