`print.anovamvl1lm` <-
function(x, ...)
    {
    fp <- format.pval(x$p.value, digits = 4)
    
    cat("\nComparisons between multivariate linear models\n")
    
    if (!is.list(x$models)) {
        cat("\nModel: ",deparse(x$models), "\n\n", sep = "")
        } else {
        cat("\nFull model:       ",deparse(x$models[[1]]), sep = "")
        cat("\nRestricted model: ",deparse(x$models[[2]]), "\n\n", sep = "")
        }
    
    cat(x$method)
    cat("\n")
    cat(paste("Q.2 = ", format(round(x$statistic,4))," with ", x$parameter, " df, p.value ", if (substr(fp, 1L, 1L) == 
        "<") fp else paste("=", fp), sep=""))
    cat("\n\n")
    invisible(x)

    }
