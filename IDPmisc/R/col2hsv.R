### col2hsv.R

col2hsv <- function(col){
  ## Author: Rene Locher
  ## Version: 2005-01-18
  rgb2hsv(col2rgb(col))
}
