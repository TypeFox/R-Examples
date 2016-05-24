########################################################################################################################
## Surcharge de la fonction print pour les objects de classe S4 VSLCMresultsContinuous et VSLCMresultsCategorical
########################################################################################################################

## Surcharge pour VSLCMresultsContinuous
setMethod(
  f="print",
  signature = c("VSLCMresultsContinuous"),
  definition = function(x){
    summary(x)    
    if (x@criteria@degeneracyrate != 1){
      cat("\n Parameters per class:\n")
      for (k in 1:x@model@g){
        if (k>1){
          cat("*******************\n")
        }
        cat("Class",k,"\n")
        cat("Proportion:",x@param@pi[k],"\n")
        tmp <- data.frame(mean=x@param@mu[,k], sd=x@param@sd[,k])
        print(tmp)
        cat("\n")
      }
    }
  }
)