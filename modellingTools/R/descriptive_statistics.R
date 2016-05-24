
#' View and return the frequency distribution of a variable.
#'
#' For continuous variables, the user can optionally specify to discretize the
#' variable into a fixed number of equal width bins, or into custom bins of the
#' user's choice. This is useful for larger datasets with many unique
#' observed values
#'
#' @param dat a tbl
#' @param var character string giving the name of the desired variable, or
#'             a single number giving the position of the desired variable
#' @param bins if 0, no discretization is performed. If a positive integer then
#'        \code{var} is binned into \code{bins} equal width ranges, and the
#'        frequency distribution of those ranges is computed. If a length > 1
#'        numeric vector, then \code{var} is binned into ranges with cutpoints
#'        defined by the unique entries of \code{bins}
#' @return a tbl containing 3 columns: \code{level} gives the unique values or
#'         bins, \code{count} gives the count in each level of \code{level} and
#'         \code{percent} gives the percentage of total observations in each
#'         level. \code{proc_freq} also automatically sends the frequency
#'         distribution to the viewer, using \code{utils::View}
#' @family descriptive
#' @details
#' R has many one-line solutions to getting the frequency distribution of a
#' variable; this function provides a unified approach that makes use of the
#' efficient data types and computation provided by the \code{dplyr} package,
#' and as a bonus, makes it easy to explore the distribution of a continuous
#' variable with many unique observations by automating discretization. The name
#' is intended to make the function more portable for SAS users who are not
#' comfortable outside their native habitat.
#'
#' @examples
#' proc_freq(faithful,"eruptions")
#' proc_freq(faithful,"eruptions",bins = 4)
#' proc_freq(faithful,"eruptions",bins = c(1,2,3,4,5))
#'
#'
#' @export
#' @import magrittr
#' @importFrom dplyr n

proc_freq <- function(dat,var,bins = 0) {

  if (!is.numeric(column_vector(dat,var))) {
    bins = 0
  } else if (length(bins) == 1 &&
             bins >= length(unique(column_vector(dat,var)))) {
    bins = 0
  }


  if (length(bins) > 1) {
    dat %<>%
      dplyr::select_(var) %>%
      dplyr::mutate(level = cut(dat[[var]],bins,include.lowest = TRUE))
  }   else if (bins > 0) {
    dat %<>%
      dplyr::select_(var) %>%
      dplyr::mutate(level = ggplot2::cut_interval(dat[[var]],n = bins))
  } else {
    dat %<>% dplyr::select_(var) %>% dplyr::mutate_(level = var)
  }

  dat %<>%
    dplyr::group_by_("level") %>%
    dplyr::summarize(count = n()) %>%
    dplyr::mutate(percent = scales::percent(count / sum(count)))

  return(dat)

}

#' Get the correlation of variables in a dataset with a given response,
#' sorted highest to lowest
#'
#' This function computes the correlation of each input variable in a dataframe
#' with a given response variable and returns a dataframe listing the variables
#' sorted in order of most to least correlated. NAs are removed from correlation
#' computations, and only numeric variables are considered.
#'
#' @param dat a tbl
#' @param response_var character string containing the name of a variable in
#'        \code{dat} that you would like the correlations to be computed with,
#'        or an integer specifying the position of this variable
#' @param parallel logical. If \code{TRUE}, parallel
#'        \code{\link[foreach]{foreach}} is used for computing correlations
#'        (if FALSE, single threaded foreach is used; still highly efficient).
#'        Default is \code{FALSE}.
#' @return a tbl with two columns: \code{var_name} gives the name of each
#'         variable and \code{correlation} gives its correlation with
#'         \code{response_var}.
#'
#' @family descriptive
#' @details
#' Use this technique for filtering out variables in the initial stages of data
#' analysis, to get more familiar with how the individual input variables relate
#' to the response variable of interest. Not recommended as a formal variable
#' selection technique, since it will ignore interactions between inputs.
#'
#' @examples
#' x <- iris
#' get_top_corrs(x,"Petal.Length")
#'
#' @export
#' @import magrittr

get_top_corrs <- function(dat,response_var,parallel = FALSE) {

  response <- column_vector(dat,response_var)
  if (!is.numeric(response)) stop("Response must be a numeric variable" +
                                    " for this to make sense. Have you" +
                                    " tried modellingTools::cross_freq?")

  dat <- dat[ ,sapply(dat,is.numeric)]

  if (parallel) {
    corrs <- foreach::foreach(nm = colnames(dat),
                              .combine = dplyr::bind_rows,
                              .multicombine = TRUE,
                              .final = function(x) {
                                stats::setNames(x,c("var_name","correlation"))
                              },
                              .packages = c("dplyr","modellingTools","stats"),
                              .export = c("nm","var_name")
    ) %dopar% {
      dplyr::data_frame(var_name = nm,
                        correlation = stats::cor(modellingTools::column_vector(dat,
                                                        nm),
                                          response,
                                          use = "na.or.complete")
      )
    }
  } else {
    corrs <- foreach::foreach(nm = colnames(dat),
                              .combine = dplyr::bind_rows,
                              .multicombine = TRUE,
                              .final = function(x) {
                                stats::setNames(x,c("var_name","correlation"))
                              },
                              .export = c("response","dat","nm","var_name")
                              ) %do% {
      dplyr::data_frame(var_name = nm,
                        correlation = stats::cor(column_vector(dat,nm),
                                          response,
                                          use = "na.or.complete")
      )
    }
  }

  corrs %<>%
    dplyr::filter(var_name != response_var) %>%
    dplyr::arrange(dplyr::desc(abs(correlation)))

  return(corrs)

}

########
# Next idea: cross_freq using foreach on the levels of a categorical variable
########


