#' Returns an unpivoted table of variable values (factor levels) associated with each branch.
#'
#' @param object an rpart object
#' @export
#' @examples
#' library(rpart)
#' fit<-rpart(Reliability~.,data=car.test.frame)
#' rpart.subrules.table(fit)
rpart.subrules.table<-function(object)  
{
  lists<-rpart.lists(object)
  leftCompares<-lapply(lists$L,function (x) attr(x,"compare"))
  rightCompares<-lapply(lists$R,function (x) attr(x,"compare"))
  leftRules<-lapply(seq_along(lists$L),function (i) setNames(data.frame(paste('L',i,sep=''),names(lists$L)[i],as.character(unlist(lists$L[i],use.names=FALSE)),NA,NA),c("Subrule","Variable","Value","Less","Greater")))
  rightRules<-lapply(seq_along(lists$R),function (i) setNames(data.frame(paste('R',i,sep=''),names(lists$R)[i],as.character(unlist(lists$R[i]),use.names=FALSE),NA,NA),c("Subrule","Variable","Value","Less","Greater")))
  
  reassign.columns<-function(object,compare)
  {
    if(grepl("<",compare))
      object$Less<-object$Value
    if(grepl(">",compare))
      object$Greater<-object$Value
    if(!grepl("=",compare))
      object$Value=NA
    return(object)
  }
  
  leftTable<-Reduce(rbind,Map(reassign.columns, leftRules, leftCompares))
  rightTable<-Reduce(rbind,Map(reassign.columns, rightRules, rightCompares))
  
  
  return(rbind(leftTable,rightTable))
}