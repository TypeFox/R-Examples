### add alpha.transparency to given color

addAlphTrsp2col <- function(col, alpha=255){
   do.call(rgb,as.list(c(col2rgb(col)[,1],max=255,alpha=alpha)))
}
