print.carbayes <- function(x,...)
{
#### Check for missingness
    if(is.null(x$samples$Y))
    {
    n.miss <- 0    
    }else
    {
    n.miss <- ncol(x$samples$Y)
    }
    
    if(class(x$localised.structure)=="list")
    {
        #### Print out the model fitted
        cat("\n#################\n")
        cat("#### Model fitted\n")
        cat("#################\n")
        cat(x$model)
        cat("Regression equation - ")
        print(x$formula)
        cat("Number of missing observations - ")
        cat(n.miss)
        cat("\n")
        
        #### Print out the results
        cat("\n############\n")
        cat("#### Results\n")
        cat("############\n")
        cat("Posterior quantities and DIC\n\n")
        print(x$summary.results)
        cat("\nDIC = ", x$modelfit[1], "     ", "p.d = ", x$modelfit[2], "     ", "LMPL = ", x$modelfit[5],"\n")

        cat("\nThe number of stepchanges identified in the random effect surface\n")
        temp <- x$localised.structure[[1]][!is.na(x$localised.structure[[1]])]
        tab <- array(NA, c(1,2))
        tab[1, ] <- c(sum(temp)/2, length(temp)/2- sum(temp)/2)
        colnames(tab) <- c("no stepchange", "stepchange")
        print(tab)
    }else if(class(x$localised.structure)=="numeric")
    {
        #### Print out the model fitted
        cat("\n#################\n")
        cat("#### Model fitted\n")
        cat("#################\n")
        cat(x$model)
        cat("Regression equation - ")
        print(x$formula)
        cat("Number of missing observations - ")
        cat(n.miss)
        cat("\n")
        
        #### Print out the results
        cat("\n############\n")
        cat("#### Results\n")
        cat("############\n")
        cat("Posterior quantities and DIC\n\n")
        print(x$summary.results)
        cat("\nDIC = ", x$modelfit[1], "     ", "p.d = ", x$modelfit[2], "     ", "LMPL = ", x$modelfit[5],"\n")
        cat("\nNumber of clusters with the number of data points in each one\n")
        print(table(paste("group", x$localised.structure, sep="")))
        
    }else
    {
        #### Print out the model fitted
        cat("\n#################\n")
        cat("#### Model fitted\n")
        cat("#################\n")
        cat(x$model)
        cat("Regression equation - ")
        print(x$formula)
        cat("Number of missing observations - ")
        cat(n.miss)
        cat("\n")
        
        #### Print out the results
        cat("\n############\n")
        cat("#### Results\n")
        cat("############\n")
        cat("Posterior quantities and DIC\n\n")
        print(x$summary.results)
        cat("\nDIC = ", x$modelfit[1], "     ", "p.d = ", x$modelfit[2], "     ", "LMPL = ", x$modelfit[5],"\n")
     }
        
return(invisible(x))        
}



