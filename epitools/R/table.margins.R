"table.margins" <-
function(x){
  #x = matrix or array
  if (!is.array(x)) 
    stop("x is not an array")  
  dd <- dim(x)
  y <- matrix(x, nrow = dd[1])
  z <- rbind(y, apply(y, 2, sum))
  y2 <- matrix(t(z), nrow = dd[2])
  z2 <- rbind(y2, apply(y2, 2, sum))
  z3 <- t(matrix(t(z2), nrow = prod(dd[-c(1, 2)])))
  fin <- array(z3, c(dd[1:2] + 1, dd[-c(1, 2)]))
  rownames(fin) <- c(rownames(x, do.NULL=FALSE),"Total")
  colnames(fin) <- c(colnames(x, do.NULL=FALSE),"Total")
  if(!is.null(names(dimnames(x)))){
    names(dimnames(fin)) <- names(dimnames(x))
  }
  fin
}
