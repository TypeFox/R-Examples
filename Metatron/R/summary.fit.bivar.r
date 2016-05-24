summary.fit.bivar<-function(object,level=0.95,...){
  x <- object
  if (!is.element("fit.bivar", class(x))) 
        stop("Argument 'x' must be an object of class \"fit.bivar\".")
   ci.level <- level
    if (ci.level <= 0 || ci.level >= 1) 
        stop("'ci.level' must be within 0 and 1")
    coeff <- x$fixef
    lvm<-x$levelmod
    vcov <- as.matrix(x$vcov)
    coeff <- as.numeric(coeff)
    coef.se <- sqrt(diag(vcov))
    tval <- coeff/coef.se
    tvalci <- qt((1 - ci.level)/2, df=x$studynum-2 , lower.tail = FALSE)
    pvalue <- 2 * (1 - pt(abs(tval),df=x$studynum-2))
    ci.lb <- coeff - tvalci * coef.se
    ci.ub <- coeff + tvalci * coef.se
    tabfixed <- cbind(coeff, coef.se, x$studynum-2, tval, pvalue, 1-ci.level, ci.lb, ci.ub)

colnames(tabfixed)<-c("coefficient","standard error", "df", "tval", "p", "alpha", "lower", "upper")

 if (is.null(x$mods)){
      cat("Fit bivariate model of Reitsma et al.(2005) without covariate.\n")
      names(x$fixef)[1:2]<-c("logit_sensitivity","logit_specificity")
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


rownames(tabfixed)<-names(x$fixef)

tabvcov<-x$vcov


rownames(tabvcov)<-colnames(tabvcov)<-names(x$fixef)
   cat("Fixed-effects coefficients:", "\n", sep = "")
   print(tabfixed)    
   cat("\n")
   cat("\n")

  cat("Covariance Matrix of Fixed Effect Parameter Estimates: \n")
  
   print(tabvcov)
   cat("\n")
}
