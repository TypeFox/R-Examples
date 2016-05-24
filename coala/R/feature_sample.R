sample_class <- R6Class("sample", inherit = feature_class,
  private = list(size = NA, ploidy = NA),
  public = list(
    initialize = function(sizes, ploidy, time) {
      assert_that(is.numeric(sizes))
      assert_that(length(sizes) >= 1)
      private$size <- sizes

      assert_that(is.number(ploidy))
      private$ploidy <- ploidy

      private$time <- private$add_parameter(time)
    },
    get_sizes = function() private$size,
    get_ploidy = function() private$ploidy,
    print = function() {
      pop <- seq(along = private$size)[private$size > 0]
      samples <- private$size[private$size > 0]
      cat("Sampling of ")
      for (i in seq(along = samples)) {
        cat(samples[i], " (pop ", pop[i], ")", sep = "")
        if (i < length(samples) - 1) cat(", ")
        else if (i == length(samples) - 1) cat(" and ")
      }
      cat(" individuals with ploidy", self$get_ploidy(),
          "at time", print_par(self$get_time()), "\n")
    }
  )
)


#' Creates a feature that represents the sampling from one population
#'
#' @param individuals The number of individuals that are sampled per population,
#'   given as a numeric vector.
#' @param ploidy The number of chromosomes that will be simulated per
#'   individual.
#' @param time The time at which the sample is taken.
#' @return The feature, which can be added to a model using `+`.
#' @keywords internal
feat_sample <- function(individuals, ploidy = 1, time = "0") {
  if (time != "0")
    stop("Samples at time different from 0 at not supported at the moment")
  sample_class$new(individuals, ploidy, time)
}


is_feat_sample <- function(feat) any("sample" == class(feat))


#' @describeIn get_features Returns a vector containing the number of
#'   haploids sampled per population. This is the default ordering of
#'   individuals used by coala.
#' @param for_sim If true, the sample size used internally for the simulation
#'   will be reported rather than the number of actual samples. The numbers
#'   can be unequal for the simulation of unphased data.
#' @export
get_sample_size <- function(model, for_sim = FALSE) {
  assert_that(is.logical(for_sim))
  sample_size <- read_cache(model, paste0("sample_size_", for_sim))

  if (is.null(sample_size)) {
    feat_mask <- sapply(model$feature, is_feat_sample)
    if (sum(feat_mask) > 1) stop("Only one sample is currently supported")
    sample_size <- model$feature[feat_mask][[1]]$get_sizes()

    if (for_sim) {
      sample_size <- sample_size * get_ploidy(model)
    } else {
      sample_size <- sample_size * get_samples_per_ind(model)
    }
    cache(model, paste0("sample_size_", for_sim), sample_size)
  }

  sample_size
}


get_ploidy <- function(model) {
  mask <- vapply(model$features, is_feat_sample, logical(1))
  if (!any(mask)) return(1L)
  if (sum(mask) > 1) stop("multiple sample features are not supported")
  feat <- model$features[mask][[1]]
  as.integer(feat$get_ploidy())
}


#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.sample <- function(feature, model) {
  sample_size <- get_sample_size(model, TRUE)
  if (length(feature$get_sizes()) == 1) return("")
  paste0("-I ", length(sample_size), " ",
         paste(sample_size, collapse = " "), " ")
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.sample <- conv_to_ms_arg.sample

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.sample <- conv_to_ms_arg.sample

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.sample <- function(feature, model) ""
