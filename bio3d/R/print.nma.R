"print.nma" <-
  function(x, nmodes=6, ...) {

    cn <- class(x)
    
    cat("\nCall:\n  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    
    cat("Class:\n  ", cn[1], " (",  cn[2], ")",
        "\n\n", sep = "")
    
    cat("Number of modes:\n  ", length(x$L), " (", x$triv.modes, " trivial)",
        "\n\n", sep="")

    inds <- (x$triv.modes+1):(x$triv.modes+nmodes)
    if(!is.null(x$frequencies)) {
      freqs <- round(x$frequencies[inds], 3)
      cat("Frequencies:\n", sep="")
    }
    else {
      freqs <- round(x$force.constants[inds], 3)
      cat("Force constants:\n", sep="")
    }
    
    i <- x$triv.modes
    for ( f in freqs ) {
      i <- i+1
      cat("  Mode ", i, ": \t", f,  "\n", sep="")
    }
    
    cat("\n")

  i <- paste( attributes(x)$names, collapse=", ")
  cat(strwrap(paste(" + attr:",i,"\n"),width=60, exdent=8), sep="\n")

    invisible(x)
  }
