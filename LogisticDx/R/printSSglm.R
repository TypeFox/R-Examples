#' @method print ss.glm
#' @export
print.ss.glm <- function(x, ...){
    d1 <- as.data.frame(x$ss)
     names(d1) <- "ss"
     print(d1)
     d2 <- as.data.frame(x$epc)
     names(d2) <- "epc"
     print(d2)
 }
