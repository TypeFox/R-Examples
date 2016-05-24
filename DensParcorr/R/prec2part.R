prec2part <- function(Precision)
{
  Precision = as.matrix(Precision)
  D.Prec = diag(diag(Precision)^(-.5))
  return(diag(2,dim(Precision)[1])-D.Prec%*%Precision%*%D.Prec)
}