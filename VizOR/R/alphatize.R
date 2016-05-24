## TODO: Make this function robust to vectors of all kinds, NAs and NULLs.

##' Apply alpha transparency to a color vector
##' 
##' @param color A vector of colors
##' @param alpha A vector of alpha values, recycled as necessary
##' @return A vector of colors, with alpha applied
##' @author David C. Norris
##' @keywords color
##' @export alphatize
alphatize <- Vectorize(function(color, alpha){
  do.call(rgb, as.list(c(col2rgb(color)[c('red','green','blue'),]/255, alpha=alpha)))
})
