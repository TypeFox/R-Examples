#' Carries out the selection step of fuzzyforest algorithm.
#'
#' Carries out the selection step of fuzzyforest algorithm.  Returns
#' data.frame with variable importances and top rated features.
#' @export
#' @param X                 A data.frame.
#'                          Each column corresponds to a feature vectors.
#'                          Could include additional covariates not a part of
#'                          the original modules.
#' @param y                 Response vector.
#' @param drop_fraction     A number between 0 and 1.  Percentage of features
#'                          dropped at each iteration.
#' @param number_selected   Number of features selected by fuzzyforest.
#' @param mtry_factor       In the case of regression, \code{mtry} is set to
#'                          \code{ceiling}(\eqn{\sqrt(p)}*\code{mtry_factor}).
#'                          In the case of classification, \code{mtry} is set to
#'                          \code{ceiling}((p/3)*\code{mtry_factor}).  If either
#'                          of these numbers is greater than p, \code{mtry} is
#'                          set to p.
#' @param min_ntree         Minimum number of trees grown in each random forest.
#' @param ntree_factor      A number greater than 1.  \code{ntree} for each
#'                          random is \code{ntree_factor} times the number
#'                          of features.  For each random forest, \code{ntree}
#'                          is set to \code{max}(\code{min_ntree},
#'                          \code{ntree_factor}*\code{p}).
#' @param num_processors    Number of processors used to fit random forests.
#' @param nodesize          Minimum nodesize
#' @return A data.frame with the top ranked features.
#' @note This work was partially funded by NSF IIS 1251151 and AMFAR 8721SC.
select_RF <- function(X, y, drop_fraction, number_selected, mtry_factor,
                      ntree_factor, min_ntree, num_processors, nodesize) {
     num_features <- ncol(X)
     selection_list <- list()
     CLASSIFICATION <- is.factor(y)
     #TUNING PARAMETER mtry_factor
     if(CLASSIFICATION == TRUE) {
       mtry <- min(ceiling(mtry_factor*sqrt(num_features)), num_features)
       if(missing(nodesize)){
         nodesize <- 1
       }
     }
     if(CLASSIFICATION == FALSE) {
       mtry <- min(ceiling(mtry_factor*num_features/3), num_features)
       if(missing(nodesize)){
         nodesize <- 5
       }
     }
     num_features <- ncol(X)
     ntree <- max(num_features*ntree_factor, min_ntree)
     target <- number_selected
     current_X <- X
     i <- 1
     while (num_features >= target){
       if(num_processors > 1) {
         rf = foreach(ntree = rep(ntree/num_processors, num_processors),
                      .combine = combine, .packages = 'randomForest') %dorng% {
                        randomForest(current_X , y, ntree = ntree, mtry = mtry,
                                     importance = TRUE, scale = FALSE, nodesize=nodesize) }
       }
       if(num_processors == 1) {
         rf <- randomForest(current_X, y, ntree = ntree, mtry = mtry,
                            importance = TRUE, scale = FALSE,
                            nodesize = nodesize)
       }
       var_importance <- rf$importance
       var_importance <- var_importance[order(var_importance[, 1],
                                              decreasing=TRUE), ]
       selection_list[[i]] <- data.frame(row.names(var_importance),
                                         round(var_importance[, 1], 4),
                                         stringsAsFactors=FALSE)
       names(selection_list[[i]]) <- c("feature_name", "variable_importance")
       i <- i + 1
       reduction <- ceiling(num_features*drop_fraction)
       if(num_features - reduction > target) {
         trimmed_varlist <- var_importance[1:(num_features - reduction), ]
         features <- row.names(trimmed_varlist)
         current_X <- current_X[, which(names(current_X) %in% features), drop=FALSE]
         num_features <- length(features)
         if(CLASSIFICATION==TRUE) {
           mtry <- min(ceiling(mtry_factor*sqrt(num_features)), dim(current_X)[2])
         }
         if(CLASSIFICATION==FALSE) {
           mtry <- min(ceiling(mtry_factor*num_features/3), dim(current_X)[2])
         }
         ntree <- max(num_features*ntree_factor, min_ntree)
       }
       else {
         num_features <- target - 1
         mod_varlist <- var_importance[, 1][1:target]
         features <- row.names(var_importance)[1:target]
         feature_list <- cbind(features, mod_varlist)
         selection_list[[i]] <- as.data.frame(cbind(features, round(mod_varlist, 4)),
                                              stringsAsFactors=FALSE)
         names(selection_list[[i]]) <- c("feature_name", "variable_importance")
       }
     }
     out <- list(feature_list, selection_list)
     return(out)
  }


#' Fits iterative random forest algorithm.
#'
#' Fits iterative random forest algorithm.  Returns
#' data.frame with variable importances and top rated features.
#' For now this is an internal function that I've used to explore how
#' recursive feature elimination works in simulations.  It may be exported at
#' a later time.
#' @param X                 A data.frame.
#'                          Each column corresponds to a feature vectors.
#' @param y                 Response vector.
#' @param drop_fraction     A number between 0 and 1.  Percentage of features
#'                          dropped at each iteration.
#' @param keep_fraction     A number between 0 and 1. Proportion features
#'                          from each module to retain at screening step.
#' @param mtry_factor       A positive number.  Mtry for each random forest
#'                          is set to
#'                          \code{ceiling}(\eqn{\sqrt{p}}\code{mtry_factor})
#'                          where \code{p} is the current number of features.
#' @param ntree_factor      A number greater than 1.  \code{ntree} for each
#'                          random is \code{ntree_factor} times the number
#'                          of features.  For each random forest, \code{ntree}
#'                          is set to \code{max}(\code{min_ntree},
#'                          \code{ntree_factor}*\code{p}).
#' @param min_ntree         Minimum number of trees grown in each random forest.
#' @param num_processors    Number of processors used to fit random forests.
#' @param nodesize          Minimum nodesize.
#' @return A data.frame with the top ranked features.
#' @note This work was partially funded by NSF IIS 1251151 and AMFAR 8721SC.
iterative_RF <- function(X, y, drop_fraction, keep_fraction, mtry_factor,
                         ntree_factor = 10, min_ntree=5000,
                         num_processors = 1, nodesize) {
  CLASSIFICATION <- is.factor(y)
  num_features <- ncol(X)
  #TUNING PARAMETER mtry_factor
  if(CLASSIFICATION == TRUE) {
    mtry <- min(ceiling(mtry_factor*sqrt(num_features)), num_features)
    if(missing(nodesize)){
      nodesize <- 1
    }
  }
  if(CLASSIFICATION == FALSE) {
    mtry <- min(ceiling(mtry_factor*num_features/3), num_features)
    if(missing(nodesize)){
      nodesize <- 5
    }
  }
  ntree <- max(num_features*ntree_factor, min_ntree)
  target <- ceiling(num_features * keep_fraction)
  current_X <- X
  while (num_features >= target){
    rf = foreach(ntree = rep(ntree/num_processors, num_processors),
                 .combine = combine, .packages = 'randomForest') %dorng% {
                 randomForest(X , y, ntree = ntree, mtry = mtry,
                 importance = TRUE, scale = FALSE, nodesize=nodesize) }
    var_importance <- importance(rf, type=1)
    var_importance <- var_importance[order(var_importance[, 1],
                                           decreasing=TRUE), ,drop=FALSE]
    reduction <- ceiling(num_features*drop_fraction)
    if(num_features - reduction > target) {
      trimmed_varlist <- var_importance[1:(num_features - reduction), ,drop=FALSE]
      features <- row.names(trimmed_varlist)
      current_X <- current_X[, which(names(current_X) %in% features), drop=FALSE]
      if(CLASSIFICATION == TRUE) {
        mtry <- min(ceiling(mtry_factor*sqrt(num_features)), num_features)
      }
      if(CLASSIFICATION == FALSE) {
        mtry <- min(ceiling(mtry_factor*num_features/3), num_features)
      }
      ntree <- max(num_features*ntree_factor, min_ntree)
    }
    else {
      num_features <- target - 1
      mod_varlist <- var_importance[, 1][1:target]
      features <- row.names(var_importance)[1:target]
      out <- cbind(features, mod_varlist)
    }
  }
  return(out)
}



