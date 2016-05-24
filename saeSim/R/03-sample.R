#' Sampling functions
#' 
#' These functions are intended to be used with \code{\link{sim_sample}} and not interactively. They are wrappers around \link[dplyr]{sample_frac} and \link[dplyr]{sample_n}.
#' 
#' @param groupVars character with names of variables to be used for grouping.
#' 
#' @inheritParams dplyr::sample_frac
#' 
#' @details \code{sample_numbers} is a vectorized version of \code{sample_number}.
#' 
#' \code{sample_cluster_number} and \code{sample_cluster_fraction} will sample clusters (all units in a cluster).
#' 
#' @rdname sampling
#' @export
#' 
#' @examples
#' sim_base_lm() %>% sim_sample(sample_number(5))
#' sim_base_lm() %>% sim_sample(sample_fraction(0.5))
#' sim_base_lm() %>% sim_sample(sample_cluster_number(5, groupVars = "idD"))
#' sim_base_lm() %>% sim_sample(sample_cluster_fraction(0.5, groupVars = "idD"))
sample_fraction <- function(size, replace = FALSE, weight = NULL, groupVars = NULL) {
  force(size); force(replace); force(weight); force(groupVars)
  function(dat) {
    if(is.null(groupVars)) {
      dat %>% sample_frac(size = size, replace = replace, weight = weight)
    } else {
      attributesToKeep <- attributes(dat)[!(names(attributes(dat)) %in% names(attributes(data.frame())))]
      dat <- group_by_(dat, groupVars) %>% 
        sample_frac(size = size, replace = replace, weight = weight) %>% as.data.frame
      attributes(dat) <- c(attributes(dat), attributesToKeep)
      dat
    }
  }
}

#' @rdname sampling
#' @export
sample_number <- function(size, replace = FALSE, weight = NULL, groupVars = NULL) {
  force(size); force(replace); force(weight); force(groupVars)
  function(dat) {
    if(is.null(groupVars)) {
      dat %>% sample_n(size = size, replace = replace, weight = weight)
    } else {
      attributesToKeep <- attributes(dat)[!(names(attributes(dat)) %in% names(attributes(data.frame())))]
      dat <- group_by_(dat, groupVars) %>% 
        sample_n(size = size, replace = replace, weight = weight) %>% as.data.frame
      
      attributes(dat) <- c(attributes(dat), attributesToKeep)
      dat
    }
  }
}

#' @rdname sampling
#' @export
sample_numbers <- function(size, replace = FALSE, groupVars = NULL) {
  mapply_by(groupVars, lapply(size, function(s) sample_number(s, replace)))
}

#' @rdname sampling
#' @export
sample_cluster_number <- function(size, replace = FALSE, weight = NULL, groupVars) {
  force(size); force(replace); force(weight); force(groupVars)
  sample_fun <- function(dat) {
    selectedGroups <- dat[groupVars] %>% unique %>% 
      sample_n(size = size, replace = replace, weight = weight) %>%
      as.data.frame
    left_join(selectedGroups, dat, by = groupVars)
  }
  preserve_attributes(sample_fun)
}

#' @rdname sampling
#' @export
sample_cluster_fraction <- function(size, replace = FALSE, weight = NULL, groupVars) {
  force(size); force(replace); force(weight); force(groupVars)
  sample_fun <- function(dat) {
    selectedGroups <- dat[groupVars] %>% unique %>% 
      sample_frac(size = size, replace = replace, weight = weight) %>%
      as.data.frame
    left_join(selectedGroups, dat, by = groupVars)
  }
  preserve_attributes(sample_fun)
}