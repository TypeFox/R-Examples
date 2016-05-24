print.dir.cond <- function(x, digits=max(3,getOption("digits") - 3),...)
  {
    cat("\nConditional Dirichlet Model MCMC output\n")
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("Parameters of the inverse gamma/normal prior:\n")
    print.default(format(x$prior,digits=digits),print.gap=2,
                  quote=FALSE)
    cat("Dirichlet precision parameter M = ",x$M,"\n\n")
    cat("Number of cycles = ", x$ncycles,"\n\n")
    if (x$start.user == TRUE) {
      cat("Starting values were by user:\n\n")
      print.default(format(x$start,digits=digits),print.gap=2,
                    quote=FALSE,...)
      }
    else if (x$start.user == FALSE) {
      nn <- length(x$start)
      cat("Default starting values were used:\n\n")
      cat("for study effects---the study estimates\n")
      cat("for mu---the average study effect = ",x$start[nn-1],"\n")
      cat("for tau---the study sd =",x$start[nn],"\n")
      }
    invisible(x)
  }

print.dir.ord <- function(x, digits=max(3,getOption("digits") - 3),...)
  {
    cat("\nOrdinary Dirichlet Model MCMC output\n")
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("Parameters of the inverse gamma/normal prior:\n")
    print.default(format(x$prior,digits=digits),print.gap=2,
                  quote=FALSE)
    cat("Dirichlet precision parameter M = ",x$M,"\n\n")
    cat("Number of cycles = ", x$ncycles,"\n\n")
    if (x$start.user == TRUE) {
      cat("Starting values were by user:\n\n")
      print.default(format(x$start,digits=digits),print.gap=2,
                    quote=FALSE,...)
      }
    else if (x$start.user == FALSE) {
      nn <- length(x$start)
      cat("Default starting values were used:\n\n")
      cat("for study effects---the study estimates\n")
      cat("for mu---the average study effect = ",x$start[nn-1],"\n")
      cat("for tau---the study sd =",x$start[nn],"\n")
      }
    invisible(x)
  }
