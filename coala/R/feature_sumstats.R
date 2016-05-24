segsites_feat_class <- R6Class("seg_sites_feat", inherit = feature_class,
  public = list(
    print = function() cat("Generating Seg. Sites\n")
  )
)

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.seg_sites_feat <- function(feature, model) {
  if (!any(vapply(get_features(model), is_feat_mutation, logical(1)))) {
    stop("model requires mutation to calculate summary statistics")
  }
  ""
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.seg_sites_feat <- conv_to_ms_arg.seg_sites_feat

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.seg_sites_feat <- conv_to_ms_arg.seg_sites_feat

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.seg_sites_feat <- conv_to_ms_arg.seg_sites_feat



trees_feat_class <- R6Class("trees_feat", inherit = feature_class,
  public = list(
    print = function() cat("Generating Trees\n")
  )
)

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.trees_feat <- function(feature, model) "-T "

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.trees_feat <- conv_to_ms_arg.trees_feat

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.trees_feat <- conv_to_ms_arg.trees_feat

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.trees_feat <- function(feature, model) {
  stop("generation of trees is not supported.")
}



files_feat_class <- R6Class("files_feat", inherit = feature_class,
  public = list(
    print = function() cat("Generating Files\n")
  )
)

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.files_feat <- ignore_par

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.files_feat <- ignore_par

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.files_feat <- ignore_par

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.files_feat <- ignore_par
