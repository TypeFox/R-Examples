
## Misc. settings and definitions

#' @import ggplot2
NULL

#' @importFrom stats cmdscale dist na.omit predict

#' @importFrom methods setClass setGeneric setMethod extends new prototype signature slot
NULL

#' @importFrom utils globalVariables
NULL

range01 <- function(x){(x-min(x))/(max(x)-min(x))} 

globalVariables(c("value", "Var1", "Var2", "X1", "X2", "object.lofscores", "MeanDecreaseAccuracy", "features", "variable", "controlobject123"))