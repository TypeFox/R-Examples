"print.geostas" <- function(x, ...) {
    
  cn <- class(x)
  ndoms <- length(unique(x$grps))
  dims <- dim(x$amsm)
  
  cat("\nCall:\n  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  
  cat("Class:\n  ", cn, "\n\n", sep = "")

  cat("Dimensions of AMSM:\n  ", dims[1], "x", dims[2], "\n\n", sep="")
  
  cat("Number of domains:\n  ", ndoms, "\n\n", sep="")

  cat("Domain size:\n")
  for(i in 1:ndoms) {
    cat("  #", i, ": ", length(x$inds[[i]]$atom), " atoms \n", sep="")
  }

  cat("\n")
  
  i <- paste( attributes(x)$names, collapse=", ")
  cat(strwrap(paste(" + attr:",i,"\n"),width=60, exdent=8), sep="\n")

  invisible(x)
}
