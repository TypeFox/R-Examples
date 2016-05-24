make.fdt_cat.multiple <- function(x,
                                  sort,
                                  decreasing)
{
  stopifnot(is.data.frame(x))

  res <- list()

  f <- sapply(x, 
              is.factor)

  for(i in 1:ncol(x)){
    if(f[i]){
      m  <- x[,i]

      fdt <- make.fdt_cat.simple(m,
                                 sort,
                                 decreasing)
      res <- c(res,
               list(fdt))
    }
  }
  valCol <- f[f]

  names(res) <- names(valCol)

  return(res)  
}
