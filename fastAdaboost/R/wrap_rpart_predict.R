#'wraps the prediction from a tree for ease of use 
#'with Rcpp
#'@noRd
#'@import rpart
#'@import stats
#'@param tree_object object of class rpart
#'@param newdata dataframe
#'@param classnames_map named vector mapping classnames to 0/1
#'@return integer_class integer vector. This contains
#'        the labels of the predicted class 
#'        converted to integers.
#'@keywords internal.
wrap_rpart_predict <- function(tree_object, newdata, classnames_map)
{
  test_learn <- predict(tree_object, newdata, type="class")
  integer_class <- ifelse( test_learn == classnames_map["A"],0, 1)
  return(integer_class)
}

#'wraps the prediction of probability from a fitted tree
#'for use with Rcpp
#'@noRd
#'@import rpart
#'@import stats
#'@param tree_object object of class rpart
#'@param newdata dataframe
#'@return pred_vec vector of probability that
#'                 example belongs to class 0
#'@keywords internal.

wrap_rpart_predict_real <-function(tree_object, newdata)
{
  pred_vec <- predict(tree_object, newdata, type="prob")[,1]
  return(pred_vec)
}