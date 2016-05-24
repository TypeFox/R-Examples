median.fdt <- function(x, ...)
{
  fdt <- with(x,
              table)

  n <- fdt[nrow(fdt), 5]

  posM <- grep(TRUE,
               n / 2 <= fdt[, 5])[1]

  brk <- with(x,
              seq(breaks['start'],
                  breaks['end'],
                  breaks['h']))

  liM <- brk[posM]

  ## It is important if 'posM' is inside of the first class
  if (posM - 1 <= 0)
    sfaM <- 0
  else
    sfaM <- fdt[(posM - 1), 5]

  fM <- fdt[posM, 2]

  h <- as.vector(with(x,
                      breaks['h']))

  res <- liM + (((n / 2) - sfaM) * h) / fM

  return(res)
}                        

