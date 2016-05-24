#'wraps the rpart call for ease of use with Rcpp
#'@noRd
#'@import rpart
#'@import stats
#'@param formula_obj formula object
#'@param newdata dataframe with requisite columns
#'@param weight_vec vector of numeric values,weights for each example
#'@param classnames_map named vector mapping classnames to 0/1.
#'@keywords internal

wrap_rpart <- function(formula_obj,newdata, weight_vec,classnames_map)
{
  formula <- as.formula(formula_obj)
  rpart_control <- rpart::rpart.control(cp=0)
  environment(formula)<-environment() #re-sets the formula environment   
  #otherwise the weight_vec is not interpreted properly
  tree_fit <- rpart::rpart(formula,newdata,weights = weight_vec, control = rpart_control)
  train_learn <- predict(tree_fit,type="class")
  train_prob <- predict(tree_fit, type="prob")[,1]
  #convert the class to 0/1 using the classname_map provided as input
  integer_class <- ifelse( train_learn == classnames_map["A"],0, 1)
  
  return_list <- list(tree = tree_fit, pred = integer_class, prob = train_prob)
  return(return_list)

}