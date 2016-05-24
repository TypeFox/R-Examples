#' Returns a list of strings summarizing the branch path to each node.
#'
#'
#' @param object an rpart object
#' @export
#' @examples
#' library(rpart)
#' fit<-rpart(Reliability~.,data=car.test.frame)
#' rpart.rules(fit)
rpart.rules<-function(object)
{
  frame<-object$frame
  ruleNums<-as.numeric(row.names(frame))  ##Convert the row names into a list of rule numbers
  is.leaf <- (frame$var == "<leaf>")
  frame[!is.leaf,"order"]<-seq_along(which(!is.leaf)) ##Number the branches to number them for matching with subrule sets
  rules<-replicate(max(ruleNums),NULL)
  rules[1]<-"NULL"
  
  ##The rule numbering convention contains the information to determine branch lineage. 
  ##Most of the potential rule numbers don't actually exist, but this will result in the creation of a NULL rule.
  for (i in as.numeric(row.names(frame))[-1])
  {
    if(i%%2==0)
    {
      rules[i]<-paste(rules[i/2],paste('L',frame[as.character(i/2),"order"],sep=''),sep=',')
    }
    else
    {
        rules[i]<-paste(rules[(i-1)/2],paste('R',frame[as.character((i-1)/2),"order"],sep=''),sep=',')
    }
  }
  rules<-lapply(rules,function (x) gsub("NULL,",'',x))
  return(rules)
}

