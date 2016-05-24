sd.fdt.multiple <- function(x, ...)
{
  res <- lapply(x,
                sd.fdt)

  return(res)
}
