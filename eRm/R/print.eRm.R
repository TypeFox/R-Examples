`print.eRm` <-
function(x,...)  {                                         #print method for all models
  cat("\n")
  cat("Results of", x$model, "estimation: \n")
  cat("\n")
  cat("Call: ", deparse(x$call), "\n")
  cat("\n")
  cat("Conditional log-likelihood:", x$loglik, "\n")
  cat("Number of iterations:", x$iter, "\n")
  cat("Number of parameters:", x$npar, "\n")
  cat("\n")
  if (x$model %in% c("RM","RSM","PCM"))                    #eta parameters
    if(is.null(x$call$W)){                                 # labelling based on whether W was specified mm 2012-05-02
      cat("Item (Category) Difficulty Parameters (eta):")  # new labelling rh 25-03-2010
    } else {
      cat("Item (Category) Parameters (eta):\nBased on design matrix W =", deparse(x$call$W))
    }
  else                                                     # now difficulty for RM, RSM, PCM
      cat("Basic Parameters eta:")
  cat("\n")
  etapar <- x$etapar
  #nameeta <- paste("eta",1:dim(x$W)[2])
  se <- x$se.eta
  result <- rbind(etapar, se)
  #colnames(result) <- nameeta
  rownames(result) <- c("Estimate", "Std.Err")
  print(result)
  cat("\n\n")
  invisible(result)
}

