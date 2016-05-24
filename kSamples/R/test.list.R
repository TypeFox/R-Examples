test.list <- function(x){
sel <- unlist(lapply(x,length))
x <- x[sel > 1]
if(length(x) <= 1) stop("need more than one block with at least 2 samples")
x
}
