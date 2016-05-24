simexit <- function(entry, all.bhr, eta.ij, x.i, max.time, pme){
  result <- NULL
  while(is.null(result)){
    u <- runif(1)
    try(result <- simexitinner(u = u, entry = entry, all.bhr = all.bhr, 
                               eta.ij = eta.ij, x.i = x.i, max.time = max.time, pme = pme), 
        silent = TRUE)
    if(is.null(result)){
      stop(paste("It was not possible to simulate a new exit time!\n", 
                 " + Check if eta functions are vectorised, and/or\n", 
                 " + increase (only those larger 0) the baseline hazard rates ", 
                 paste(names(all.bhr), collapse = ", "), ", and/or (if used)\n", 
                 " + decrease partial-markov effects!", sep = ""))
    }
  }
  return(result)
}