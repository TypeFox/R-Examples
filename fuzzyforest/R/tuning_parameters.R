#' Set Parameters for Screening Step of
#' Fuzzy Forests
#'
#' Creates \code{screen_control} object for
#' controlling how feature selection
#' will be carried out on each module.
#' @export
#' @param drop_fraction     A number between 0 and 1.  Percentage of features
#'                          dropped at each iteration.
#' @param keep_fraction     A number between 0 and 1. Proportion of features
#'                          from each module that are retained from screening step.
#' @param mtry_factor       In the case of regression, \code{mtry} is set to
#'                          \code{ceiling}(\eqn{\sqrt(p)}*\code{mtry_factor}).
#'                          In the case of classification, \code{mtry} is set to
#'                          \code{ceiling}((p/3)*\code{mtry_factor}).  If either
#'                          of these numbers is greater than p, \code{mtry} is
#'                          set to p.
#'
#' @param ntree_factor      A number greater than 1.  \code{ntree} for each
#'                          random forest is \code{ntree_factor} times the number
#'                          of features.  For each random forest, \code{ntree}
#'                          is set to \code{max}(\code{min_ntree},
#'                          \code{ntree_factor}*p).
#' @param min_ntree         Minimum number of trees grown in each random forest.
#' @return An object of type screen_control.
#' @references
#' Daniel Conn, Tuck Ngun, Christina M. Ramirez (2015). Fuzzy Forests: a New
#' WGCNA Based Random Forest Algorithm for Correlated, High-Dimensional Data,
#' Journal of Statistical Software, Manuscript in progress.
#' @examples
#' drop_fraction <- .25
#' keep_fraction <- .1
#' mtry_factor <- 1
#' min_ntree <- 5000
#' ntree_factor <- 5
#' screen_params <- screen_control(drop_fraction=drop_fraction,
#'                                 keep_fraction=keep_fraction,
#'                                 mtry_factor=mtry_factor,
#'                                 min_ntree=min_ntree,
#'                                 ntree_factor=ntree_factor)
#' @note This work was partially funded by NSF IIS 1251151 and AMFAR 8721SC.
screen_control <- function(drop_fraction=.25, keep_fraction=.05,
                           mtry_factor=1, min_ntree=5000, ntree_factor=10) {
  obj <- list()
  obj$drop_fraction <- drop_fraction
  obj$keep_fraction <- keep_fraction
  obj$mtry_factor <- mtry_factor
  obj$min_ntree <- min_ntree
  obj$ntree_factor <- ntree_factor
  class(obj) <- "screen_control"
  return(obj)
}

#' Set Parameters for Selection Step of
#' Fuzzy Forests
#'
#' Creates \code{selection_control} object for
#' controlling how feature selection
#' will be carried out after features from different
#' modules have been combined.
#' @export
#' @param drop_fraction     A number between 0 and 1.  Percentage of features
#'                          dropped at each iteration.
#' @param number_selected   A positive number. Number of features
#'                          that will be selected by fuzzyforests.
#' @param mtry_factor       In the case of regression, \code{mtry} is set to
#'                          \code{ceiling}(\eqn{\sqrt(p)}*\code{mtry_factor}).
#'                          In the case of classification, \code{mtry} is set to
#'                          \code{ceiling}((p/3)*\code{mtry_factor}).  If either
#'                          of these numbers is greater than p, \code{mtry} is
#'                          set to p.
#' @param ntree_factor      A number greater than 1.  \code{ntree} for each
#'                          random forest is \code{ntree_factor} times the number
#'                          of features.  For each random forest, \code{ntree}
#'                          is set to \code{max}(\code{min_ntree},
#'                          \code{ntree_factor}*\code{p}).
#' @param min_ntree         Minimum number of trees grown in each random forest.
#' @return An object of type selection_control.
#' @references
#' Daniel Conn, Tuck Ngun, Christina M. Ramirez (2015). Fuzzy Forests: a New
#' WGCNA Based Random Forest Algorithm for Correlated, High-Dimensional Data,
#' Journal of Statistical Software, Manuscript in progress.
#' @examples
#' drop_fraction <- .25
#' number_selected <- 10
#' mtry_factor <- 1
#' min_ntree <- 5000
#' ntree_factor <- 5
#' select_params <- select_control(drop_fraction=drop_fraction,
#'                                 number_selected=number_selected,
#'                                 mtry_factor=mtry_factor,
#'                                 min_ntree=min_ntree,
#'                                 ntree_factor=ntree_factor)
#' @note This work was partially funded by NSF IIS 1251151 and AMFAR 8721SC.
select_control <- function(drop_fraction=.25, number_selected=5,
                              mtry_factor=1, min_ntree=5000, ntree_factor=10) {
  obj <- list()
  obj$drop_fraction <- drop_fraction
  obj$number_selected <- number_selected
  obj$mtry_factor <- mtry_factor
  obj$min_ntree <- min_ntree
  obj$ntree_factor <- ntree_factor
  class(obj) <- "select_control"
  return(obj)
}


#' Set Parameters for WGCNA Step of
#' Fuzzy Forests
#'
#' Creates \code{WGCNA_control} object for
#' controlling WGCNA will be carried out.
#' @export
#' @param power             Power of adjacency function.
#' @param ...               Additional arguments.
#'                          See \code{\link[WGCNA]{blockwiseModules}} for
#'                          details.
#' @return An object of type WGCNA_control.
#' @references
#' Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene
#' Co-Expression Network Analysis", Statistical Applications in Genetics and
#' Molecular Biology: Vol. 4: No. 1, Article 17
#' @examples
#' WGCNA_params <- WGCNA_control(p=7, minModuleSize=30, TOMType = "unsigned",
#'                               reassignThreshold = 0, mergeCutHeight = 0.25,
#'                               numericLabels = TRUE, pamRespectsDendro = FALSE)
#' @note This work was partially funded by NSF IIS 1251151.
WGCNA_control <- function(power=6, ...) {
  extra_args <- list(...)
  obj <- list()
  obj$power <- power
  obj$extra_args <- extra_args
  class(obj) <- "WGCNA_control"
  return(obj)
}







