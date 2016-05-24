simexitinner <- function(u, entry, all.bhr, eta.ij, x.i, max.time, pme){
  suppressWarnings(new.exit <- uniroot(unirootf, interval = c(entry, max.time * 1e4),
                      u = u, entry = entry, all.bhr = all.bhr, eta.ij = eta.ij, x.i = x.i, pme = pme)$root)
  return(list(new.exit = new.exit, u = u))}