 
 #############################
 #### print.summary.ldblm ####
 #############################

print.summary.ldblm <-function(x,digits = 2,...){
 
 if (!inherits(x,"summary.ldblm")) 
   stop("use only with \"summary.ldblm\" objects")
 
 # print the call 
 x$call[[1]]<-as.name("ldblm")
 cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")

 # print the family
 if (!is.character(x$family))
   cat(gettextf("\nfamily: %s",x$family$family),"\n") 
 else
   cat("\nfamily: gaussian\n")

 # print the Residuals
 cat("\nResiduals:\n")
 print(summary(as.numeric(format(round(x$residuals,2)))))

 # print the Sum of Squares
 cat(gettextf("\nR-squared : %s",format(round(x$r.squared,digits=4))))
 
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

