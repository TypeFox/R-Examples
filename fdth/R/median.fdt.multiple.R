median.fdt.multiple <- function(x, ...)
{
  res <- lapply(x,
                median.fdt)

  return(res)
}
