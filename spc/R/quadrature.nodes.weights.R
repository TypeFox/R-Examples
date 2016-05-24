quadrature.nodes.weights <- function(n, type="GL", x1=-1, x2=1) {
  if ( n < 1 )			stop("n has to be a natural number")
  qtyp <- pmatch(type, c("GL", "Ra")) - 1
  if ( is.na(qtyp) )		stop("invalid quadrature type")
  if ( x1 >= x2 )		stop("x1 must be smaller than x2")
  
  nw <- .C("quadrature_nodes_weights",
           as.integer(n),
           as.double(x1), as.double(x2),
           as.integer(qtyp),
           ans=double(length=2*n), PACKAGE="spc")$ans
           
  qnw <- data.frame(nodes=nw[1:n], weights=nw[-(1:n)])         
  qnw
}