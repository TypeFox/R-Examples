rtnlc <-
function (m, x) {
  # function for returning listcolumn x in a matrix of lists 
  # ("ReTurN List Column")
  # m = inMatrix
  aFun = function(r, c, m, x) {as.numeric(unlist(strsplit(m[r, c],",")))[x]}
  vect_aFun = Vectorize(aFun, vectorize.args = c('r','c'))
  return(outer(1:nrow(m), 1:ncol(m) , FUN = vect_aFun, m, x))
}
