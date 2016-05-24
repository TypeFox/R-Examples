summary.anova<-function (object,ndec=NULL,...) {
 if (is.null(ndec)) ndec=5
 if (class(object)=="anova.RPm") {
    cat("     - SUMMARY anova.RPm - \n")
    cat("\n p-value for Bonferroni method \n" )
    print(round(object$p.Bonf,ndec))
    cat("\n  p-value for False Discovery Rate method \n")
     print(round(object$p.FDR,ndec))
    if (!is.null(object$p.Boot)){
           cat("\n p-value for Bootstrap method \n" )
           print(round(object$p.Boot,ndec))
     }
 }
 if (class(object)=="anova.hetero") {
     cat("\n  - SUMMARY ANOVA HETEROCEDASTHIC - \n")
          print(object)
  }
cat("\n")
output<-object
}


