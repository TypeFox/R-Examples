### #######
### General
### #######

divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")

## ######
## Binary
## ######

print.binIRT <- function(x, ...) {
    cat("\n")
    cat(divider)
    cat("\n")
    cat(paste("Call:\n",
              paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
              sep = ""
              )
        )

    cat(paste("\tIndividuals: ", x$n, "\n",
              "\tItems: ", x$j, "\n",
              "\tLatent Dimensions: ", x$d, "\n",
              sep = ""
              )
        )

    cat("\n")

    cat(paste("\tIterations: ", x$runtime$iters, "\n",
              "\tTolerance: ", x$runtime$tolerance, "\n",
              "\tConvergence Status: ", as.integer(x$runtime$conv),"\n",
              "\tClassification Success Rate: ", round(x$fit$csr, 4), "\n",
              "\tGMP: ", round(x$fit$gmp, 4), "\n",
              sep = ""
              )
        )
    cat("\n")
    cat(divider)
    cat("\n")

    return(invisible(x))
}

### #######
### Ordinal
### #######

print.ordIRT <- function(x, ...) {
    cat("\n")
    cat(divider)
    cat("\n")
    cat(paste("Call:\n",
              paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
              sep = ""
              )
        )

    cat(paste("\tIndividuals: ", x$N, "\n",
              "\tItems: ", x$J, "\n",
              "\tLatent Dimensions: ", 1, "\n",
              sep = ""
              )
        )

    cat("\n")

    cat(paste("\tIterations: ", x$runtime$iters, "\n",
              "\tTolerance: ", x$runtime$tolerance, "\n",
              "\tConvergence Status: ", x$runtime$conv ,"\n",
              sep = ""
              )
        )
    cat("\n")
    cat(divider)
    cat("\n")

    return(invisible(x))
}

### #######
### Dynamic
### #######

print.dynIRT <- function(x, ...) {
    cat("\n")
    cat(divider)
    cat("\n")
    cat(paste("Call:\n",
    paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
    sep = ""
    )
    )
    
    cat(paste("\tIndividuals: ", x$N, "\n",
    "\tItems: ", x$J, "\n",
    "\tLatent Dimensions: ", 1, "\n",
    sep = ""
    )
    )
    
    cat("\n")
    
    cat(paste("\tIterations: ", x$runtime$iters, "\n",
    "\tTolerance: ", x$runtime$tolerance, "\n",
    "\tConvergence Status: ", x$runtime$conv ,"\n",
    sep = ""
    )
    )
    cat("\n")
    cat(divider)
    cat("\n")
    
    return(invisible(x))
}

### #######
### Hierarchical
### #######

print.hierIRT <- function(x, ...) {
    cat("\n")
    cat(divider)
    cat("\n")
    cat(paste("Call:\n",
    paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
    sep = ""
    )
    )
    
    cat(paste("\tGroups (G): ", x$G, "\n",
    "\tMembers (I): ", x$N$I, "\n",
    "\tItems (J): ", x$N$J, "\n",
    "\tTotal observed votes (L): ", x$N$L, "\n",
    "\tLatent Dimensions (D): ", x$N$D, "\n",
    sep = ""
    )
    )
    
    cat("\n")
    
    cat(paste("\tIterations: ", x$runtime$iters, "\n",
    "\tTolerance: ", x$runtime$tolerance, "\n",
    "\tConvergence Status: ", x$runtime$conv ,"\n",
    sep = ""
    )
    )
    cat("\n")
    cat(divider)
    cat("\n")
    
    return(invisible(x))
}


