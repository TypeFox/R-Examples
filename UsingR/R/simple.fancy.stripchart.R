##' function to make a fancier stripchart
##'
##' @param l data
##' @return NULL
##'
##' @export
"simple.fancy.stripchart" <-
  function(l) {
    stripchart(l,pch=1,group.names=names(l))
    n = length(l)
    for(i in 1:n) {
      abline(i,0,lty=3)
      points(mean(l[[i]]),i-.5/n,pch=17,cex=2)
    }
  }
