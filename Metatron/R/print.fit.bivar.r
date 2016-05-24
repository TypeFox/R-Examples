print.fit.bivar<-function (x,...) 
{
   if (!is.element("fit.bivar", class(x))) 
        stop("Argument 'x' must be an object of class \"fit.bivar\".")

   names(x$fixef)[1:2]<-c("logit_sensitivity","logit_specificity") 

   lvm<-x$levelmod 
  aleatorio<-matrix(c(x$rancoef[1,1],x$rancoef[2,1],x$rancoef[1,2],x$rancoef[2,2]),nrow=2)


  if (is.null(x$mods)){
      cat("Fit bivariate model of Reitsma et al.(2005) without covariate.\n")
      
}
  else if (!is.null(x$bi.both)){
      cat("Fit bivariate model of Reitsma et al.(2005) with covariate",x$mf.mods,"that affects both sensitivity and specificity.\n")   
    for (i in 1:length(lvm)){
 names(x$fixef)[i]<-paste0("logit_sensitivity_",lvm[i])
   names(x$fixef)[length(lvm)+i]<-paste0("logit_specificity_",lvm[i])
  }

}
  else if (!is.null(x$bi.sen)){
      cat("Fit bivariate model of Reitsma et al.(2005) with covariate",x$mf.mods,"that affects only sensitivity.\n")
   names(x$fixef)[1]<-paste0("logit_specificity_",lvm[1])
   names(x$fixef)[2]<-paste0("logit_sensitivity_",lvm[1])   
   for (i in 3:(length(lvm)+1)){
 names(x$fixef)[i]<-paste0("logit_sensitivity_",lvm[i-1])
     } 
}
  else if (!is.null(x$bi.spe)){
      cat("Fit bivariate model of Reitsma et al.(2005) with covariate",x$mf.mods,"that affects only specificity.\n")
   names(x$fixef)[1]<-paste0("logit_sensitivity_",lvm[1])
   names(x$fixef)[2]<-paste0("logit_specificity_",lvm[1])
   for (i in 3:(length(lvm)+1)){
 names(x$fixef)[i]<-paste0("logit_specificity_",lvm[i-1])
     }
}
 cat("\n")
 cat("Fixed-effects coefficients:", "\n", sep = "")
 table <- formatC(x$fixef, format = "f")
 table <- gsub("NA", " -", table)
    print(table, quote = FALSE, right = TRUE, print.gap = 2)
    cat("\n")
    cat("\n")
 cat("Random-effects coefficients which model the between-study variability,shown as variance-covariance matrix:", "\n", sep = "")
  
       rownames(aleatorio)<-colnames(aleatorio)<-c("logit_sensitivity","logit_specificity")   

 print(aleatorio)
    cat("\n")
    cat("\n")
    cat(x$studynum, " studies ", x$obsnum ," classifications ", x$df.fixef, " fixed and 3 random-effects parameters", "\n", sep = "")
    table <- c(x$logLik, x$AIC, x$BIC)
    names(table) <- c("logLik", "AIC", "BIC")
    table <- formatC(table, format = "f")
    print(table, quote = FALSE, right = TRUE, print.gap = 2)
    cat("\n")   
}

