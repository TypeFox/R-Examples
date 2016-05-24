#'unpivot data matrix
#'@export
#'@param m matrix or dataframe
#'@examples
#'x = matrix(1:25,ncol=5)
#'x = as.data.frame(x)
#'colnames(x) = letters[1:5]
#'rownames(x) = LETTERS[1:5]
#'unpivot(x)
unpivot = function(m){
  m = as.matrix(m)
  m = data.frame( rep(rownames(m), ncol(m)),
                  rep(colnames(m), each=nrow(m)),
                  c(m) )
  return(m)
}
