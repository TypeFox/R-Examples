#' Returns an unpivoted table of branch paths (subrules) associated with each node.
#'
#'
#' @param object an rpart object
#' @export
#' @examples
#' library(rpart)
#' fit<-rpart(Reliability~.,data=car.test.frame)
#' rpart.rules.table(fit)
rpart.rules.table<-function(object)
{
  rules<-rpart.rules(object)
  ff<-object$frame
  ff$rules<-unlist(rules[as.numeric(row.names(ff))])
  ruleList<-lapply(row.names(ff),function (name) setNames(data.frame(name,
                                                                     (strsplit(ff[name,'rules'],split=',')),
                                                                     ff[name,'var']=="<leaf>"
                                                                     ),
                                                          c("Rule","Subrule","Leaf")))
  combinedRules<-Reduce(rbind,ruleList)
  
  return(combinedRules)
  
}