"print.enma" <-
  function(x, ...) {

    cn      <- class(x)
    nstruct <- nrow(x$fluctuations)
    dims    <- dim(x$U.subspace)

    if(is.null(x$call$rm.gaps))
      rm.gaps <- TRUE
    else if(x$call$rm.gaps=="T" || x$call$rm.gaps=="TRUE")
      rm.gaps <- TRUE
    else
      rm.gaps <- FALSE
    
    if(is.null(x$call$fit))
      fit <- TRUE
    else if(x$call$fit=="T" || x$call$fit=="TRUE")
      fit <- TRUE
    else
      fit <- FALSE


    
    cat("\nCall:\n  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    
    cat("Class:\n  ", cn,
        "\n\n", sep = "")
    
    cat("Number of structures:\n  ", nstruct,
        "\n\n", sep="")

    cat("Attributes stored:\n")
    
    if(!is.null(x$full.nma))
      cat("  - Full 'nma' objects\n")
    if(!is.null(x$rmsip))
      cat("  - Root mean square inner product (RMSIP)\n")
    if(!is.null(x$fluctuations))
      cat("  - Aligned atomic fluctuations\n")
    
    if(rm.gaps)
      cat("  - Aligned eigenvectors (gaps removed)\n")
    else
      cat("  - Aligned eigenvectors (gaps not removed)\n")
    
    cat("  - Dimensions of x$U.subspace: ", dims[1L], "x", dims[2L], "x", dims[3L], sep="")

    cat("\n\n")

    if(fit)
      cat("Coordinates were aligned prior to NMA calculations")
    else
      cat("Coordinates were NOT aligned prior to NMA calculations")
    
    cat("\n\n")

  i <- paste( attributes(x)$names, collapse=", ")
  cat(strwrap(paste(" + attr:",i,"\n"),width=60, exdent=8), sep="\n")

    invisible(x)
  }
