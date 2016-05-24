#'Print adaboost.m1 model summary
#'
#'S3 method to print an adaboost object
#'
#'Displays basic information on the model,
#'such as function call, dependent variable,
#'the number of trees, and weights assigned to 
#'each tree
#'
#'@param x object of class adaboost
#'@param ... arguments passed to print.default
#'@return None
#'@export
#'@examples
#'fakedata <- data.frame( X=c(rnorm(100,0,1),rnorm(100,1,1)), Y=c(rep(0,100),rep(1,100) ) )
#'fakedata$Y <- factor(fakedata$Y)
#'test_adaboost <- adaboost(Y~X, fakedata, 10)
#'print(test_adaboost)
#'
#'@seealso \code{\link{print.real_adaboost}}

print.adaboost <-function(x,...)
{
  object <- x
  print(object$call)
  print(object$formula)
  cat("Dependent Variable: ", object$dependent_variable,"\n",sep="")
  cat("No of trees:",length(object$trees),"\n",sep="")
  cat("The weights of the trees are:",object$weights,"\n",sep="")

}


#'Print real adaboost model summary
#'
#'S3 method to print a real_adaboost object
#'
#'Displays basic information on the model,
#'such as function call, dependent variable
#'and the number of trees
#'
#'@param x object of class real_adaboost
#'@param ... arguments passed to print.default
#'@return None
#'@export
#'@examples
#'fakedata <- data.frame( X=c(rnorm(100,0,1),rnorm(100,1,1)), Y=c(rep(0,100),rep(1,100) ) )
#'fakedata$Y <- factor(fakedata$Y)
#'test_real_adaboost<- real_adaboost(Y~X, fakedata, 10)
#'print(test_real_adaboost)
#'@seealso \code{\link{print.adaboost}}

print.real_adaboost <-function(x,...)
{
  object <- x
  print(object$call)
  print(object$formula)
  cat("Dependent Variable: ", object$dependent_variable,"\n",sep="")
  cat("No of trees:",length(object$trees),"\n",sep="")
}