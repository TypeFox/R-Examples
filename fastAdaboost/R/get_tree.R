#'Fetches a decision tree
#'
#'returns a single weak decision tree classifier which is part
#'of the strong classifier
#'
#'returns an individual tree from the adaboost object
#'This can provide the user with some clarity on the 
#'individual building blocks of the strong classifier
#'
#'@param object object of class adaboost
#'@param tree_num integer describing the tree to get
#'@return object of class rpart
#'@export
#'@examples 
#'fakedata <- data.frame( X=c(rnorm(100,0,1),rnorm(100,1,1)), Y=c(rep(0,100),rep(1,100) ) )
#'fakedata$Y <- factor(fakedata$Y)
#'test_adaboost <- adaboost(Y~X, fakedata, 10)
#'tree <- get_tree(test_adaboost,5)
#'@seealso \code{\link{adaboost}}

get_tree <-function(object, tree_num)
{
  if(class(object)[1]!="adaboost")
    stop("object is not of class adaboost",call.=F)
  tree_list <- object$trees
  if(tree_num> length(tree_list))
  {
    stop(paste("requested tree number",tree_num,
               "is greater than no. of adaboost iterations:",length(tree_list)), call.=F)
  }
  return(tree_list[tree_num])
}