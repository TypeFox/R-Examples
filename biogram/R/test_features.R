#' Permutation test for feature selection
#'
#' Performs a feature selection on positioned n-gram data using a Fisher's 
#' permutation test.
#'
#' @param target \code{integer} vector with target information (e.g. class labels).
#' @param features \code{integer} matrix of features with number of rows equal 
#' to the length of the target vector.
#' @param criterion criterion used in permutation test. See \code{\link{calc_criterion}} 
#' for the list of possible criterions.
#' @param adjust name of p-value adjustment method. See \code{\link[stats]{p.adjust}}
#' for the list of possible values. If \code{NULL}, no adjustment is done.
#' @param threshold \code{integer}. Features that occur less than \code{threshold}
#' and more often than \code{nrow(features)-threshold} are discarded from the permutation 
#' test.
#' @param quick \code{logical}, if \code{TRUE} Quick Permutation Test (QuiPT) is used.
#' @param times number of times procedure should be repeated. Ignored if \code{quick} is 
#' \code{TRUE}.
#' 
#' @details Since the procedure involves multiple testing, it is advisable to use one
#' of the avaible p-value adjustment methods. Such methods can be used directly by 
#' specifying the \code{adjust} parameter.
#' @return an object of class \code{\link{feature_test}}.
#' @note Both \code{target} and \code{features} must be binary, i.e. contain only 0 
#' and 1 values.
#' 
#' Features occuring too often and too rarely are considered not informative and may be 
#' removed using the threshold parameter.
#' @export
#' @keywords nonparametric
#' @references 
#' Radivojac P, Obradovic Z, Dunker AK, Vucetic S, 
#' \emph{Feature selection filters based on the permutation test} in 
#' Machine Learning: ECML 2004, 15th European 
#' Conference on Machine Learning, Springer, 2004.
#' @seealso 
#' \code{\link{calc_criterion}} - computes selected criterion.
#' 
#' \code{\link{distr_crit}} - distribution of criterion used in QuiPT.
#' 
#' \code{\link{summary.feature_test}} - summary of results.
#' 
#' \code{\link{cut.feature_test}} - aggregates test results in groups based on feature's
#' p-value.
#' @examples
#' #significant feature
#' tar_feat1 <- create_feature_target(10, 390, 0, 600) 
#' #significant feature
#' tar_feat2 <- create_feature_target(9, 391, 1, 599)
#' #insignificant feature
#' tar_feat3 <- create_feature_target(198, 202, 300, 300)
#' test_res <- test_features(tar_feat1[, 1], cbind(tar_feat1[, 2], tar_feat2[, 2], 
#'                           tar_feat3[, 2]))
#' summary(test_res)
#' cut(test_res)
test_features <- function(target, features, criterion = "ig", adjust = "BH", 
                          threshold = 1, quick = TRUE, times = 1e5) {
  
  valid_criterion <- check_criterion(criterion)
  
  #criterion function
  crit_function <- function(target, features)
    calc_criterion(target, features, valid_criterion[["crit_function"]])
  
  #few tests for data consistency
  if (!all(target %in% c(0, 1))) {
    stop("target is not {0,1}-valued vector")
  }
  if (nrow(features) != length(target)) {
    stop("target and feature have different lengths")
  }
  
  apply(features, 2, function(feature) {
    if (!all(feature %in% c(0,1)) ) {
      stop("feature is not {0,1}-valued vector")
    }
  })
  
  feature_size <- if (class(features) == "simple_triplet_matrix") {
    slam::col_sums(features)
  } else {
    colSums(features)
  }
  
  #eliminate non-infomative features
  features <- features[, feature_size > threshold & feature_size < 
                         (nrow(features) - threshold)]
  
  p_vals <- if(quick) {
    
    # compute distribution once
    feature_size <- unique(feature_size)
    
    dists <- lapply(feature_size, function(i){
      t <- create_feature_target(i, abs(sum(target) - i), 0, 
                                 abs(length(target) - sum(target))) 
      distr_crit(t[, 1], t[, 2], criterion = criterion)
    })
    
    names(dists) <- feature_size
    
    apply(features, 2, function(feature) {
      feature <- as.matrix(feature, ncol = 1)
      n <- length(target)
      estm <- crit_function(target, feature)
      dist <- dists[[paste(sum(feature))]]
      1 - dist[which.max(dist[, "criterion"] >= estm - 1e-15), "cdf"]
    })
  } else {
    #slow version
    rowMeans(crit_function(target, features) <= 
               replicate(times, crit_function(sample(target), features)))
  }
  
  #p-values sometimes are a tiny little bit bigger than 1
  p_vals[p_vals > 1] <- 1
  
  if(!is.null(adjust))
    p_vals <- p.adjust(p_vals, method = adjust)
  
  create_feature_test(p_value = p_vals, 
                      criterion = valid_criterion[["nice_name"]],
                      adjust = adjust,
                      times = ifelse(quick, NA, times),
                      occ = calc_occurences(target, features))
}

#calculates occurences of features in target+ and target- groups
calc_occurences <- function(target, features) {
  target_b <- as.bit(target)
  len_target <- length(target)
  pos_target <- sum(target)
  occ <- apply(features, 2, function(i)
    fast_crosstable(target_b, len_target, pos_target, i))[c(1, 3), ]/
    c(pos_target, len_target - pos_target)
  rownames(occ) <- c("pos", "neg")
  occ
}
