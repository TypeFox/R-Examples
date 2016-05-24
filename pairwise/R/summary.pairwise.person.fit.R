#' @export summary.pairwise.person.fit
#' @title S3 Summary for Person-Fit-Statistics
#' @description S3 summary method for object of class\code{c("pairwise.item.fit", "data.frame" )}
#' @param object object of class\code{"pairwise.person.fit", "data.frame" }
#' @param sort logical with default \code{sort=FALSE} - if set to \code{sort=TRUE} persons are ordered by absolute FIT. 
#' @param by character passing the type of Fit-Statistic to sort by - ignored when \code{sort=FALSE}. valid options are: \code{"INFIT.ZSTD"} (default), \code{"OUTFIT.MSQ"}, \code{"OUTFIT.ZSTD"} and \code{"INFIT.MSQ"}.
#' @param decreasing see \code{\link{order}}
#' @param ... other parameters passed trough - see \code{\link{order}}

########################### hier die summary method fuer pairwise.person.fit #############################
summary.pairwise.person.fit<-function(object, sort=FALSE, by="INFIT.ZSTD", decreasing=FALSE, ...){
  fit_1 <- object[ , c(1:3,8:11)]
  colnames(fit_1)[4:7] <- c("OUTFIT.MSQ","OUTFIT.ZSTD","INFIT.MSQ","INFIT.ZSTD" )
  
  if(sort==TRUE){  
    result <- fit_1[order( abs(fit_1[,by]),decreasing = decreasing),]
    cat("(ordered by absolute value of" , by, ")","\n")
    }else{result <- fit_1}
    return(result)
  }
