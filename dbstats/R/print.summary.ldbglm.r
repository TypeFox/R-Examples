 
 ##############################
 #### print.summary.ldbglm ####
 ##############################

print.summary.ldbglm <-function(x,digits = 2,...){
 
 if (!inherits(x,"summary.ldbglm")) 
   stop("use only with \"summary.ldbglm\" objects")
 
 # print the call 
 x$call[[1]]<-as.name("ldbglm")
 cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")

 # print the family
 if (!is.character(x$family))
   cat(gettextf("\nfamily: %s",x$family$family),"\n") 
 else
   cat("\nfamily: gaussian\n")
 
 # Deviance Residuals
 cat("\nDeviance Residuals:\n")
 print(summary(as.numeric(format(round(x$deviance.resid,digits=4)))))
 
 # dispersion
 if (x$family$family %in% c("poisson","binomial"))
   cat(gettextf("\n(Dispersion parameter for %s family taken to be %i)",
                x$family$family,x$dispersion),"\n")
 else
   cat(gettextf("\n(Dispersion parameter for %s family taken to be %f)",
                x$family$family,x$dispersion),"\n")
 
 # deviance: null or estimated model
 cat(gettextf("\n    Null deviance: %s  on %i degrees of freedom",
              format(round(x$null.deviance,digits)),x$df.null))
 cat(gettextf("\nResidual deviance: %s  on %s equivalent number of degrees of freedom",
              format(round(x$residual.deviance,digits=digits)),format(round(x$df.residual,digits=digits))),"\n")
 
 # print the Number of Observations
 cat(gettextf("\nNumber of Observations:    %i",x$nobs))

 # print the Trace of smoother matrix
 cat(gettextf("\nTrace of smoothing matrix: %s",format(round(x$trace.hat,2))),"\n") 
   
  cat("\nSummary of distances between data:\n")
  print(x$summary.dist1) 
 
  # print the used bandwidth 
  if (x$method.h!="user.h")
   cat(gettextf("\nOptimal bandwidth h : %f",x$h.opt),"\n") # print h.opt  
  else 
   cat(gettextf("\nUser bandwidth h : %f",x$h.opt),"\n") # print h.opt  
 
  cat(paste(gettextf("Percentile of bandwidth in the distance matrix= %s",
               format(round(x$percentile.h.opt,2))),"%",sep=""),"\n\n") 
  # print the appropriate statistic according to the using method.h 
  if(!is.null(x$crit.value)){
    cat("Bandwidth choice based on ",x$method.h,"\n")
    cat(paste(x$method.h, "value criterion :", format(round(x$crit.value,digits)),"\n"))
  }
 
 # print the using kernel 
 cat(gettextf("\nKind of kernel= %s",x$kind.kernel),"\n")  
 cat("\n")
}

