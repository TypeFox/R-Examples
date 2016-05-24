#' @rdname LLSRxy
#' @title Dataset aggregation Function
#' @description This function returns a data.frame from given two unidimensionals vectors, XC and YX.
#' @details The function returns a data.frame after merging two unidimensional vectors.
#' @export LLSRxy
#' @param YC - Component Y's concentration in the TOP Y-rich phase.
#' @param XC - Component X's concentration in the BOTTOM X-rich phase.
#' @return XYdt - Standard bidimensional data.frame used in most of functions available in this package.
#' @examples
#' #
#' Xdt<-peg4kslt[,1]
#' #
#' Ydt<-peg4kslt[,1]
#' #
#' LLSRxy(Xdt,Ydt)
LLSRxy <- function(XC,YC) {
  # convert variables into vectors
  xc <- as.vector(XC)
  yc <- as.vector(YC)
  # and combine them into a dataframe
  XYdt <- data.frame(XC = xc,YC = yc)
  #return data silently - should it be Visible or hidden?
  invisible(XYdt)
}
