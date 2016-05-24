mfv.fdt.multiple <- function(x, ...)
{
  res <- lapply(x,
                mfv.fdt)

  return(res)
}
