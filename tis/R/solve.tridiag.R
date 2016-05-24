solve.tridiag <- function(a, b, ...){
  a <- as.matrix(a)
  da <- dim(a)
  if(da[1] != da[2]) 
    stop(paste(deparse(substitute(a)), "is not a square matrix"))
  n <- as.integer(da[1])
  if(missing(b)) 
    b <- diag(n)
  ra <- row(a)
  ca <- col(a)
  storage.mode(a) <- storage.mode(b) <- "double"
  .C("dgtsv", n, 
     nrhs = dim(as.matrix(b))[2], 
     dl = a[ra == ca +1], 
     d = a[ra == ca], 
     du = a[ra == ca - 1], 
     x = b, 
     ldb = n, 
     info = as.integer(0),
     PACKAGE = "tis")$x
}
