"print.sse" <-
  function(x, ...) {

    cn <- class(x)

    cat("\n")
    cat("Call:\n  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")

    cat("Class:\n  ", cn, "\n\n", sep = "")

    cat("Helices: ", length(x$helix$start), "")
    if(length(x$helix$start)>0) {
      j <- 0
      for(i in 1:length(x$helix$start)) {
        if(j%%5==0)
          cat("\n  ")
        
        tmpout <- paste(x$helix$start[i], "-", x$helix$end[i],
                        " (", x$helix$chain[i], ")", sep="")
        cat(format(tmpout, justify="right", width=15))
        j <- j+1
      }
    }

    cat("\n\n")
    cat("Sheets:  ", length(x$sheet$start), "")
    if(length(x$sheet$start)>0) {
      j <- 0
      for(i in 1:length(x$sheet$start)) {
        if(j%%5==0)
          cat("\n  ")
        
        tmpout <- paste(x$sheet$start[i], "-", x$sheet$end[i],
                      " (", x$sheet$chain[i], ")", sep="")
        cat(format(tmpout, justify="right", width=15))
        j <- j+1
      }
    }
    
    cat("\n\n")
    cat("Turns:  ", length(x$turn$start), "")
    if(length(x$turn$start)>0) {
      j <- 0
      for(i in 1:length(x$turn$start)) {
        if(j%%5==0)
          cat("\n  ")
        
        tmpout <- paste(x$turn$start[i], "-", x$turn$end[i],
                      " (", x$turn$chain[i], ")", sep="")
        cat(format(tmpout, justify="right", width=15))
        j <- j+1
      }
    }

    cat("\n\n")

    if(is.null(x$call$resno))
      resno <- TRUE
    else if(x$call$resno=="T" || x$call$resno=="TRUE")
      resno <- TRUE
    else
      resno <- FALSE

    if(resno)
      cat("Output is provided in residue numbers")
    else
      cat("Output is provided in residue (sequential) identifiers")

    cat("\n\n")
    invisible(x)
  }
