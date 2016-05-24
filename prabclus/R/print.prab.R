"print.prab" <-
function(x, ...){
  cat("Presence-absence matrix object with\n")
  cat(x$n.species," species and ",x$n.regions," regions,\n")
  cat("including ")
  nbind <- FALSE
  for (i in 1:length(x$nb)) if (!identical(x$nb[[i]],numeric(0))) nbind <- TRUE
  if (nbind & x$nbbetweenregions) cat("regions neighborhoods and \n")
  else{
    if (nbind & !x$nbbetweenregions)
      cat("species/individuals neighborhoods and \n")
    else
      cat("\n")
  }
  cat("between-species distance matrix of type ",x$distance,".\n")
}
