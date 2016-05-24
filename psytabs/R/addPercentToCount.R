addPercentToCount <-
function(x, y) {
  pastePercentToNumeric <- function(x, y) paste(x, " (", y, ")", sep="")
  data.frame(mapply(pastePercentToNumeric, x, y))
}
