#' @export summary.pers
#' @title S3 Summary for Thetas
#' @description S3 summary method for object of class\code{"pers"}
#' @param object object of class\code{"pers"}
#' @param short logical with default \code{short=TRUE} - if set to \code{short=FALSE} a "data.frame" with WLE estimates (and their respective standard errors) for every row (person) in the original dataset will be returned.
#' @param sortwle logical wether to order persons by ability - ignored when \code{short=TRUE}
#' @param ... other parameters passed trough

########################### hier die summary method fuer pers #############################
summary.pers<-function(object, short=TRUE, sortwle=FALSE, ...){
  personen <- object$pers
  items <- object$pair
  result <- as.data.frame(personen[1:6])
  rownames(result) <- personen[["persID"]]
  
  if(short==FALSE){  
    if(sortwle==TRUE){ 
    result <- result[order(result[,5] ),]  
    cat("(ordered by Theta) \n")
    }
    cat("WLE Estimates and SE by Persons: \n") 
    #print(result)  
    return(result)
  }
  
  if(short==TRUE){
    erg <- (aggregate(result[,c(1,3:6)],by=list(result$book,result$NA.group,result$raw),FUN=unique))
    #(aggregate(result[,c(1,3:6)],by=list(result$book,result$NA.group,result$raw),FUN=table))
        cat("Person estimates by 'booklet', 'NA.group' and Scoregroup:","\n","\n")
    print(erg[,4:8])
    cat("\n", "WLE Reliability:",(object$WLE.rel$r.WLE.rel), "; ( N =",object$WLE.rel$n.WLE.rel,")", "\n","(",(dim(object$pers)[1]-object$WLE.rel$n.WLE.rel),"persons without WLE estimate )" , "\n","(",(object$WLE.rel$N.perf),"persons with perfect response vectors excluded for WLE reliability estimation )","\n")  
  }  
}
