`print.kin.cohort` <-
function(x,descriptive=TRUE, cumrisk=TRUE,hazard=FALSE,survival=FALSE, logrank=TRUE, HR=TRUE, digits=5,...){
   cat("Kin-cohort analysis\n")
   if(inherits(x,"chatterjee")){
     cat("\nMarginal Likelihood method\n")
     if(descriptive){
        cat("\nPerson Years table\n")
        print(x$events,na="")
        cat("\nAllele frequency =",x$f,"\n")
     }
     if(cumrisk){
        cat("\nCumulative Risk\n")
        print(round(x$cumrisk,digits))
     }
     if(survival){
        cat("\nSurvival\n")
        print(round(1-x$cumrisk[,1:2],digits))
     }
     if(hazard){
        cat("\nHazard\n")
        print(round(x$hazard,digits))
     }
     if(HR){
        cat("\nAverage Hazard Ratio\n")
        print(round(exp(x$logHR),digits))
     }
     if(logrank &!is.null(x$logrank)){
        cat("\nLogrank test p-value\n")
        print(signif(x$logrank, 2))
     }
   }
   if(inherits(x,"wacholder")){
     cat("\nWacholder Moments method\n")
     if(descriptive){
        cat("\nPerson Years table\n")
        print(x$events,na="")
        cat("\nAllele frequency =",x$f,"\n")
     }
     if(cumrisk){
        cat("\nCumulative Risk\n")
        print(round(x$cumrisk,digits))
     }
     if(survival){
        cat("\nSurvival\n")
        print(round(1-x$cumrisk[,1:2],digits))
     }
     if(logrank &!is.na(x$logrank)){
        cat("\nLogrank test p-value=", format(signif(x$logrank, 2)),"\n")
     }
  }
}

