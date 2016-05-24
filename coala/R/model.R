#' Create a Coalescent Model
#'
#' This creates a basic coalescent model to which more features, loci,
#' parameters and summary statistics can be added later. Data under the model
#' can be simulated using the \code{\link[=simulate.coalmodel]{simulate}} function.
#'
#' @inheritParams feat_sample
#' @param sample_size Defines the number of populations and the number of
#'   individual sampled from each population. Given as an integer vector where
#'   each entry gives the number of individuals sampled from the corresponding
#'   population.
#' @param loci_number You can optionally add a number of loci with equal length
#'   to the model. This gives to number of loci to add.
#' @param loci_length This gives the length of the loci to add.
#' @return The basic coalescent model which can be extended with features,
#'   parameters, loci and summary statistics.
#' @export
#' @seealso The 'coala-intro' vignette for a general description on how to extend
#'          models.
#' @seealso For checking which simulators can be used for this model:
#'          \code{\link{check_model}}
#' @seealso For adding mutation or for a list of other features:
#'          \code{\link{feat_mutation}}
#' @seealso For adding loci: \code{\link{locus_single}},
#'          \code{\link{locus_averaged}}, \code{\link{locus_trio}}
#' @seealso For a generating DNA sequences or for a list of summary statistics:
#'          \code{\link{sumstat_dna}}
#' @importFrom assertthat assert_that
#' @examples
#' # A model with one population and 20 unlinked loci:
#' model <- coal_model(10, 20) +
#'   feat_mutation(5) +
#'   sumstat_tajimas_d()
#' check_model(model)
#' simulate(model)
#'
#' # A model with two populations:
#' model <- coal_model(c(13, 18), 5) +
#'   feat_migration(.5, symmetric = TRUE) +
#'   sumstat_trees()
#' check_model(model)
#' simulate(model)
#'
#' # A model with 10 populations:
#' model <- coal_model(rep(2, 10), 5) +
#'   feat_migration(.5, symmetric = TRUE) +
#'   sumstat_trees()
#' check_model(model)
#' simulate(model)
#'
#' # A model with recombination:
#' model <- coal_model(20, 1, 1000) +
#'   feat_recombination(10) +
#'   feat_mutation(5) +
#'   sumstat_four_gamete()
#' check_model(model)
#' simulate(model)
coal_model <- function(sample_size, loci_number = 0,
                       loci_length = 1000, ploidy = 1) {

  model <- list(features = list(),
                loci = list(),
                parameter = list(),
                sum_stats = create_sumstat_container(),
                has_variation = FALSE,
                scaling_factor = 1,
                id = get_id())

  class(model) <- c("coalmodel", "coalmodelpart")

  # Add sample sizes
  model <- model + feat_sample(sample_size, ploidy)

  # Add locus
  if (loci_number > 0) {
    model <- model + locus_averaged(loci_number, loci_length)
  }

  model
}


is.model <- function(model) {
  "coalmodel" %in% class(model)
}

model_part <- R6Class("coalmodelpart")


# Selects a program for simulation that is capable of all current features
select_simprog <- function(model) {
  name <- read_cache(model, "simprog")

  if (is.null(name)) {
    priority <- -Inf

    for (simprog_name in ls(simulators)) {
      simprog <- get_simulator(simprog_name)
      valid <- try(simprog$get_cmd(model), silent = TRUE)
      if (all(class(valid) != "try-error")) {
        if (simprog$get_priority() > priority) {
          name <- simprog
          priority <- simprog$get_priority()
        }
      }
    }

    if (is.null(name)) warning("No suitable simulation software found!")
    cache(model, "simprog", name)
  }

  name
}


add_variation <- function(model) {
  model$has_variation <- TRUE
  model
}


has_variation <- function(model) model$has_variation


has_trios <- function(model) {
  sum(get_locus_length_matrix(model)[ , c(1:2, 4:5)]) > 0
}


get_snp_positions <- function(seg_sites, model, relative=TRUE) {
  lapply(1:length(seg_sites), function(locus) {
    pos <- get_positions(seg_sites[[locus]])
    locus_length <- get_locus_length(model, locus, total = FALSE)

    # Nothing changes without trios
    if (length(locus_length) == 1) {
      if (relative) return(pos)
      else return(pos * locus_length)
    }

    # Convert if we have trios
    trio_locus <- get_trio_locus(seg_sites[[locus]])
    if (is.null(trio_locus)) trio_locus <- 0
    pos[trio_locus == -1] <- pos[trio_locus == -1] * locus_length[1]
    pos[trio_locus == 0] <- pos[trio_locus == 0] * locus_length[3] +
      sum(locus_length[1:2])
    pos[trio_locus == 1] <- pos[trio_locus == 1] * locus_length[5] +
      sum(locus_length[1:4])
    if (relative) pos <- pos / sum(locus_length)
    pos
  })
}
