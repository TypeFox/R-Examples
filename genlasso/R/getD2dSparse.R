getD2dSparse <- function(dim1, dim2) {
  D1 = bandSparse(dim1*dim2, m=dim1*dim2, k=c(0,1), diagonals=list(rep(-1,dim1*dim2),rep(1,dim1*dim2-1)))
  D1 = D1[(Seq(1,dim1*dim2) %% dim1) != 0,]
  D2 = bandSparse(dim1*dim2-dim1, m=dim1*dim2, k=c(0,dim1), diagonals=list(rep(-1,dim1*dim2),rep(1,dim1*dim2-1)))
  return(rBind(D1,D2))
}
