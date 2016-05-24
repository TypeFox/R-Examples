# Check chosen criterion
# 
# Checks if the criterion is viable or matches it to the list of implemented 
# criterions.
# 
# @param input_criterion a \code{character} string, criterion from input.
# @param criterion_names list of implemented criterions, always in lowercase.
# @export
# @return a list of three: 
# \itemize{
# \item{criterion name,}
# \item{its function,}
# \item{nice name for outputs.}
# }
# @seealso
# Calculate the value of criterion: \code{\link{calc_criterion}}.
check_criterion <- function(input_criterion, criterion_names = c("ig", "kl")) {
  #think twice about grep
  valid_name <- criterion_names[grepl(tolower(input_criterion), criterion_names)]
  
  if (length(valid_name) == 0)
    stop("Name ", input_criterion, " cannot be associated with any available criterion.")
  
  if (length(valid_name) > 1)
    stop("Name ", input_criterion, " is too ambiguous. Rerun with more precise name.")
  
  
  criterion_data <- switch(valid_name,
                           ig = list(crit_function = calc_ig, nice_name = "Information Gain"),
                           kl = list(crit_function = calc_kl, nice_name = "Kullback-Leibler divergence"))
  #TO DO - should also return the full name of criterion for purpose of summaries/plots
  c(crit_name = valid_name, criterion_data)
}

#' Calculate criterion
#'
#' Calculates independently chosen statistical criterion for each feature versus target vector.
#'
#' @details Permutation test implemented in \code{biogram} uses several criterions to filter 
#' important features. Each can be used by \code{\link{test_features}} by specifying
#' \code{criterion} parameter.
#' 
#' Possible criterions are:
#' \describe{
#'   \item{ig}{Information Gain. Calculated using \code{\link{calc_ig}}.}
#'   \item{kl}{Kullback-Leibler divergence. Calculated using \code{\link{calc_kl}}.}
#' }
#' @inheritParams test_features
#' @param criterion_function a function calculating criterion. For a full list, see See also.
#' @return a \code{integer} vector of length equal to the number of features 
#' containing computed information gain values.
#' @note Both \code{target} and \code{features} must be binary, i.e. contain only 0 
#' and 1 values.
#' @seealso \code{\link{test_features}}.
#' @export
#' @examples 
#' calc_criterion(sample(0L:1, 100, replace = TRUE), 
#'                matrix(sample(0L:1, 400, replace = TRUE), ncol = 4),
#'                calc_ig)
calc_criterion <- function(target, features, criterion_function) {
  tar_bit <- as.bit(target)
  l_tar <- length(target)
  pos_tar <- sum(target)
  props_tar <- c(l_tar - pos_tar, pos_tar)/l_tar
  #entrophy
  ES <- - sum(props_tar * entlog(props_tar))
  apply(features, 2, function(single_feature) 
    criterion_function(single_feature, tar_bit, l_tar, pos_tar, ES))
}