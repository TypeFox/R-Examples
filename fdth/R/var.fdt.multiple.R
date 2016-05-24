var.fdt.multiple <- function(x, ...)
{
  res <- lapply(x,
                var.fdt)

  return(res)
}
