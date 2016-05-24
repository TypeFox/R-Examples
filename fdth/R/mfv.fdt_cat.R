mfv.fdt_cat <- function(x, ...)
{
  fdt <- x

  y <- fdt[, 2]

  posMFV <- grep(max(y),
                 y)

  res <- fdt[posMFV, 2]

  names(res) <- fdt[posMFV, 1]

  return(res)
}                        
