transp <- function (col, alpha = 0.5) 
{
  res <- apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255, 
                                                c[3]/255, alpha))
  return(res)
}