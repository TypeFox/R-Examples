"print.pca" <-
  function(x, nmodes=6, ...) {

    cn <- class(x)
    
    cat("\nCall:\n  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    
    cat("Class:\n  ", cn, "\n\n", sep = "")
    
    cat("Number of eigenvalues:\n  ", length(x$L), 
        "\n\n", sep="")

    inds <- 1:nmodes
    e <- round(x$L[inds], 3)
    p <- (x$L[inds]/sum(x$L)) * 100
    d <- data.frame( "Eigenvalue"=e, 
                     "Variance"=round(p,3), 
                     "Cumulative"=round(cumsum(p),3), 
                     row.names = paste("   PC",inds))

    #cat("Eigenvalues:\n", sep="")
    print(d)

    cat("\n",paste("  (Obtained from", nrow(x$z), "conformers with", ncol(x$z), "xyz input values)."))
    i <- paste( attributes(x)$names, collapse=", ")
    cat("\n",strwrap(paste(" + attr:",i,"\n"),width=60, exdent=8), sep="\n")

    invisible(d)
  }

