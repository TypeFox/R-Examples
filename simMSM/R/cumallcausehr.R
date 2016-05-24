cumallcausehr <- function(entry, exit, all.bhr, eta.ij, x.i, pme){
  result <- integrate(allcausehr, lower = entry, upper = exit, all.bhr = all.bhr,
                      eta.ij = eta.ij, x.i = x.i, pme = pme, subdivisions = 10000)
  return(result$value)}