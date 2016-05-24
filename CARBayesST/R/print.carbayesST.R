print.carbayesST <- function(x,...)
{
     if(class(x$localised.structure)=="list")
     {
     #### Print out the model fitted
     cat("\n#################\n")
     cat("#### Model fitted\n")
     cat("#################\n")
     cat(x$model)
     cat("Regression equation - ")
     print(x$formula)

     #### Print out the results
     cat("\n############\n")
     cat("#### Results\n")
     cat("############\n")
     cat("Posterior quantities for selected parameters and DIC\n\n")
     print(x$summary.results)
     cat("\nDIC = ", x$modelfit[1], "     ", "p.d = ", x$modelfit[2], "     ", "LMPL = ", x$modelfit[5], "\n")
     cat("\nThe number of stepchanges identified in the random effect surface")
     cat("\nthat satisfy Prob(w_ij < 0.5|data) > 0.99 is \n")
     temp <- x$localised.structure[[2]][!is.na(x$localised.structure[[2]])]
     tab <- array(NA, c(1,2))
     tab[1, ] <- c(sum(temp)/2, (length(temp)- sum(temp))/2)
     colnames(tab) <- c("stepchange", "no stepchange")
     print(tab)
     #print(sum(temp))
     #cat(" which is ")
     #print(round(100 * mean(temp),2))
     #cat("%\n")
     }else if(class(x$localised.structure)=="numeric")
     {
     #### Print out the model fitted
     cat("\n#################\n")
     cat("#### Model fitted\n")
     cat("#################\n")
     cat(x$model)
     cat("Regression equation - ")
     print(x$formula)

     #### Print out the results
     cat("\n############\n")
     cat("#### Results\n")
     cat("############\n")
     cat("Posterior quantities for selected parameters and DIC\n\n")
     print(x$summary.results)
     cat("\nDIC = ", x$modelfit[1], "     ", "p.d = ", x$modelfit[2], "     ", "LMPL = ", x$modelfit[5], "\n")
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

     #### Print out the results
     cat("\n############\n")
     cat("#### Results\n")
     cat("############\n")
     cat("Posterior quantities for selected parameters and DIC\n\n")
     print(x$summary.results)
     cat("\nDIC = ", x$modelfit[1], "     ", "p.d = ", x$modelfit[2], "     ", "LMPL = ", x$modelfit[5], "\n")
     }
     
return(invisible(x))
}