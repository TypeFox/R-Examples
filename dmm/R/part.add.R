part.add <-
function(x,p,prow,pcol){
# part.add() - add a partition p to a matrix x at prow,pcol
  x[prow : (prow + nrow(p) -1), pcol : (pcol + ncol(p) -1)] <- p[,]
  return(x)
}
