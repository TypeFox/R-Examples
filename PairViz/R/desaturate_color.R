desaturate_color <- function(cols,frac=.8){
 rgbcols <- col2rgb(cols,TRUE)
 cols <- rgb2hsv(rgbcols[1:3,])
 cols[2,] <- cols[2,]*frac
 cols <- hsv(cols[1,],cols[2,],cols[3,],alpha=rgbcols[4,]/255)
 return(cols)
}
