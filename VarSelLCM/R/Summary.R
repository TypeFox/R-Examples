########################################################################################################################
## Surcharge de la fonction summary pour les objects de classe S4 VSLCMresultsContinuous et VSLCMresultsCategorical
########################################################################################################################

## Surcharge pour VSLCMresultsContinuous
setMethod(
  f="summary",
  signature = c("VSLCMresultsContinuous"),
  definition = function(object){
    cat("Data set:\n   Number of individuals:", object@data@n,"\n")
    cat("   Number of continuous variables:", object@data@d, "\n")
    val <- round(100*(1-mean(object@data@notNA)),2)
    if (val>0)
      cat("   Percentile of missing values:", ,"\n\n")
    cat("Model:\n   Number of components:", object@model@g, "\n   Number of relevant variables for the clustering",sum(object@model@omega),"\n")
    if (sum(object@model@omega)>0){
      cat("\nNames of the relevant variables for the clustering:\n  ")
      print(colnames(object@data@data)[which(object@model@omega==1)])
    }   
    if (object@criteria@degeneracyrate != 1){
      cat("\n\nInformation Criteria:\n")
      cat("   loglike:", object@criteria@loglikelihood,"\n")    
      cat("   BIC:    ", object@criteria@BIC,"\n")    
      cat("   ICL:    ", object@criteria@ICL,"\n") 
      if (object@strategy@vbleSelec)
        cat("   MICL:   ", object@criteria@MICL,"\n")      
    }
    cat("\n")
    if (object@criteria@degeneracyrate>0.1)
      cat("Warnings:\n  The rate of degeneracy for the EM algorithm is", object@criteria@degeneracyrate,"\n" )
    
  }
)