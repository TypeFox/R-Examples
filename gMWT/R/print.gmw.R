`print.gmw` <- function(x, ...){
 X <- list()
 X$p.values <- x$p.values
 print(X$p.values,...)
}