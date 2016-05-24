plot.whatif <- function(x, type = "f", numcf = NULL, eps = FALSE, ...)  {

    #USER ERROR INPUT CHECKS
    #Check if plot argument is valid
    if (!((identical(type, "f") || identical(type, "l") || identical(type, "b"))))  {
      stop("argument 'type' must be one of the characters ''l'', ''f'', or ''b''")
    }

    #Check if numcf is a numeric vector without missing data if
    #user supplies the argument
    if (!is.null(numcf))  {
      if (!(is.vector(numcf) && is.numeric(numcf)))  {
        stop("argument 'numcf' must be a numeric vector")
      }
      na.fail(numcf)
    }

    #LOCAL VARIABLES
    x.val <- as.numeric(dimnames(x$cum.freq)[[2]])
   
    #LOCAL FUNCTIONS
    do.plot <- function(x, type, x.val, i)  {
      if (x$in.hull[i])  {
        if (identical(type, "f") || identical(type, "b"))  {
          lines(x.val, x$cum.freq[i,], lty = 1)
        }
        if (identical(type, "l") || identical(type, "b"))  {
          lines(lowess(x.val, x$cum.freq[i,],
            f = 0.3), lty = 1, col = 4)
        }
      }  else  {
        if (identical(type, "f") || identical(type, "b"))  {
          lines(x.val, x$cum.freq[i,], lty = 2)
        }
        if (identical(type, "l") || identical(type, "b"))  {
          lines(lowess(x.val, x$cum.freq[i,],
            f = 0.3), lty = 2, col = 4)
        }
      }
    }

    #PLOTTING
    if (eps)  {
      if (is.null(numcf))  {
        filename <- paste("graph_", type, "_all", ".eps", sep="")
      }  else {
        filename <- paste("graph_", type, "_", numcf[1], ".eps", sep="")
      }
      postscript(file=filename, onefile=FALSE, horizontal=FALSE)
    }
    plot(x.val, x$cum.freq[1,], xlab =
      "Distance", ylab =
      "Cumulative % of data within distance", ylim = c(0, 1),
      xlim = c(min(x.val), max(x.val)), type = "n", main = 
      "Cumulative Frequencies of Distances to Counterfactual")
    if (is.null(numcf))  {
      for (i in 1:dim(x$cum.freq)[1])  {
        do.plot(x, type, x.val, i)
      }
    }  else  {
      for (j in 1:length(numcf)) {
        do.plot(x, type, x.val, numcf[j])
      }
    }
    if (eps) {
      dev.off()
    }
}
