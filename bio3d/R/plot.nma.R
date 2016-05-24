"plot.nma" <-
  function(x, pch = 16, col = par("col"), cex = 0.8,
           mar = c(6, 4, 2, 2), ...) {
    
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    
    par(cex = cex, mar = mar)
    layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
    
    if(!is.null(x$frequencies)) {
      freqs <- x$frequencies
      main <- "Frequencies"
    }
    else {
      freqs <- x$force.constants
      main <- "Force constants"
    }

    if(length(freqs)>=100)
      n <- 100
    else
      n <- length(freqs)
    
    plot(x$L[1:n], type = "h", pch = pch, xlab = "Mode index",
         ylab = "", main = "Eigenvalues", col = col)

    plot(freqs[1:n], type = "h", pch = pch, xlab = "Mode index",
         ylab = "", main = main, col = col) 
    
    plot.bio3d(x$fluctuations, pch = pch, xlab = "Residue index",
         ylab = "", main = "Fluctuations", col = col, ...)
  }
