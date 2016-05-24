
#' @title Multinomial Sufficient Statistics
#' @description Calculate observed-data sufficient statistics, marginally-observed summary statistics,
#' or enumerate all possible observed patterns from a multivariate multinomial dataset. 
#' @param dat A \code{data.frame}. All variables must be factors.
#' @param output A string specifying the desired output. One of \code{c("x_y", "z_Os_y", "possible.obs")}.
#' \code{"x_y"} indicates the observed-data sufficient statistics, \code{"z_Os_y"} indicates the
#' marginally-observed summary statistics, and \code{"possible.obs"} indicates the possible observed 
#' patterns.
#' @return A \code{data.frame} containing either sufficient statistics or possible observed patterns.
#' 
#' @examples \dontrun{
#'  data(tract2221)
#'  obs_suff_stats <- multinomial_stats(tract2221, output= "x_y")
#'  marg_obs_suff_stats <- multinomial_stats(tract2221, output= "z_Os_y")
#' }
#' @export
multinomial_stats <- function(dat, output= c("x_y", "z_Os_y", "possible.obs")) {
  # 00. validate inputs
  #----------------------------------------------
  output <- match.arg(output, several.ok= FALSE)
  
  if (!all(apply(dat, 2, is.factor))) {
    # enforce factor variables
    dat <- data.frame(apply(dat, 2, function(x) as.factor(x)))
  }
  
  # 01. build empty contingency table
  #----------------------------------------------
  if (output != "z_Os_y") {
    enum_comp <- expand.grid(sapply(dat, function(x) return(levels(x))))
  } else {
    enum <- expand.grid(sapply(dat, function(x) return(c(levels(x), NA))))
    enum_comp <- enum[stats::complete.cases(enum),] 
    enum_miss <- enum[!stats::complete.cases(enum),]
    enum_miss <- enum_miss[apply(enum_miss, 1, function(x) !all(is.na(x))),] # not all missing
  }
  rownames(enum_comp) <- 1:nrow(enum_comp) # y \in Y
  
  # 02. get counts / sufficient statistics
  # return outputs
  #----------------------------------------------
  if (output == "x_y") {
    # complete data sufficient statistics
    dat_comp <- dat[stats::complete.cases(dat),]
    return(count_levels(dat_comp, enum_list= enum_comp, hasNA= "no"))
  } else if (output == "z_Os_y") {
    dat_miss <- dat[!stats::complete.cases(dat),]
    # missing data marginal sufficient statistics
    return(count_levels(dat_miss, enum_list= enum_miss, hasNA= "count.miss"))
  } else if (output == "possible.obs") {
    # just return possible combinations, no counts
    return(enum_comp)
  }
}
  