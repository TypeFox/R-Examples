quantile.fdt.multiple <- function(x, ...)
{
  res <- lapply(x,
                quantile.fdt, ...)

  return(res)
}
