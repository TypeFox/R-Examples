mean.fdt.multiple <- function(x, ...)
{
  res <- lapply(x,
                mean.fdt, ...)

  return(res)
}
