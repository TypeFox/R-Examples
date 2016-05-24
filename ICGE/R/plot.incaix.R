plot.incaix <- function(x, y, ...){
# x must be a list with three components
   k <- length(x$well_class)
   p <- (x$well_class)/(x$Ni_cluster)*100
   barplot(t(p), names.arg=1:k, ylab="% well classified units", xlab="Group")
}
