matrix2<-function(x, stages, byrow=TRUE)
{
   if(!is.vector(x)){ x<-unlist(x)}
   rows <- sqrt(length(x))
   if(ceiling(rows) != floor(rows)){stop("x will not convert to a square matrix")}
   if(missing(stages)){stages <- 1:rows}
   if(length(stages) != rows){stages <- stages[1:rows]} # or exit with error
   y<-matrix(x, nrow=rows, byrow=byrow, dimnames=list(stages, stages))
   #class(y) <- c("projection", "matrix")   # add class?
   y
}
