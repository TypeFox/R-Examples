quantile.fdt <- function(x,
                         ...,
                         i=1, 
                         probs=seq(0, 1, 0.25))
{
  fdt <- with(x,
              table)

  n <- fdt[nrow(fdt), 5]

  posQ <- grep(TRUE,
               i * n / (length(probs) - 1) <= fdt[, 5])[1]

  brk <- with(x,
              seq(breaks['start'],
                  breaks['end'],
                  breaks['h']))

  liQ <- brk[posQ]

  ## It is important if 'posQ ' is inside of the first class
  if (posQ - 1 <= 0)
    sfaQ <- 0
  else
    sfaQ <- fdt[(posQ - 1), 5]

  fQ <- fdt[posQ, 2]

  h <- as.vector(with(x,
                      breaks['h']))

  res <- liQ + (((i * n / (length(probs) - 1)) - sfaQ) * h) / fQ

  return(res)
}                        
