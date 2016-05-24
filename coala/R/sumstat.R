#' Base Class for Summary Statistics
#'
#' If you want to create additional summary statistics for coala, create
#' \pkg{R6} classes that inherit from this object.
#'
#' @include model.R
#' @keywords internal
#' @export
sumstat_class <- R6Class("sumstat", inherit = model_part,
  private = list(
    name = NA,
    req_files = FALSE,
    req_trees = FALSE,
    req_segsites = FALSE,
    transformation = NULL
  ),
  public = list(
    initialize = function(name, transformation) {
      assert_that(is.character(name))
      assert_that(length(name) == 1)
      private$name <- name

      assert_that(is.function(transformation))
      private$transformation <- transformation
    },
    calculate = function(seg_sites, trees, files, model) {
      stop("Overwrite this function with the calculation of the statistic.")
    },
    check = function(model) {
      # Optional functions that checks if a model is compatible with the stat.
      # Should throw an informative error if not.
      # This gets executed before the stat is added to the model.
      invisible(TRUE)
    },
    get_name = function() private$name,
    requires_files = function() private$req_files,
    requires_segsites = function() private$req_segsites,
    requires_trees = function() private$req_trees,
    print = function() cat(class(self)[1], "\n"),
    transform = function(x) private$transformation(x)
  )
)


is.sum_stat <- function(sumstat) inherits(sumstat, "sumstat")


# Written to `model$sum_stats` on initialization
create_sumstat_container <- function() list()


# Add a summary statistic to a model
add_to_model.sumstat <- function(sum_stat, model, feat_name) {
  sum_stat$check(model)

  if (sum_stat$get_name() %in% names(model$sum_stats))
    stop("Can not add ", feat_name, " to model: ",
         "There is already a statistic with name ", sum_stat$get_name())

  if (sum_stat$requires_files() && !requires_files(model))
    model <- model + files_feat_class$new()
  if (sum_stat$requires_segsites() && !requires_segsites(model))
    model <- model + segsites_feat_class$new()
  if (sum_stat$requires_trees() && !requires_trees(model))
    model <- model + trees_feat_class$new()

  # Save the statistic
  model$sum_stats[[sum_stat$get_name()]] <- sum_stat

  # Update cache
  model$id <- get_id()
  model
}


#' @param pop The population for which aspects are returned
#' @describeIn get_features Returns the summary statistics in the model
#' @export
get_summary_statistics <- function(model) {
  model$sum_stats
}


calc_sumstats <- function(model, segsites_list = NULL, trees = NULL,
                          files = NULL, ...) {

  assert_that(is.model(model))
  if (requires_segsites(model)) {
    assert_that(!is.null(segsites_list))
    assert_that(all(vapply(segsites_list, is_segsites, logical(1))))

    if (has_ign_singletons(model)) {
      segsites_list <- remove_singletons(segsites_list)
    }
  }
  if (requires_trees(model)) assert_that(!is.null(trees))
  if (requires_files(model)) assert_that(!is.null(files))

  stats <- lapply(model$sum_stats, function(stat) {
      stat$transform(stat$calculate(segsites_list, trees, files, model))
  })

  c(stats, list(...))
}


calc_sumstats_from_sim <- function(seg_sites, trees, files, model,
                                   pars, cmds, simulator) {

  if (missing(pars)) pars <- numeric(0)
  assert_that(is.model(model))

  if (is.list(cmds) && simulator$get_name() != "seqgen") {
    cmds <- do.call(c, cmds)
  }

  # Process seg_sites for trios and unphase if neccessary
  if (requires_segsites(model)) {
    if (has_trios(model) && simulator$get_name() != "seqgen") {
      seg_sites <- conv_for_trios(seg_sites, model)
    }

    if (is_unphased(model)) {
      seg_sites <- unphase_segsites(seg_sites,
                                    get_ploidy(model),
                                    get_samples_per_ind(model))
    }
  }

  calc_sumstats(model, seg_sites, trees, files,
                pars = pars, cmds = cmds, simulator = simulator$get_info())
}


#' Calculate summary statistics for biological data
#'
#' This function calculates a model's summary statistic from biological data.
#' The data needs to be provided as a list of segregating sites objects. These
#' objects can be create using the \code{\link{create_segsites}} function.
#'
#' @param model The coala model. The summary statistics present in this model
#'   will be calculated. The model should fit to the data, in particular
#'   regarding the number of loci and haploids.
#' @param segsites_list Either a list of \code{segsites} objects, or an object
#'   that can be converted using \code{\link{as.segsites}}. It is possible
#'   to specify additional argument for the conversion using the \code{...}
#'   argument.
#' @param tree_list Not yet implemented.
#' @param trios If your model is using locus trios, then you
#'   can create these by combining individual loci. This is a list that defines
#'   which loci are combined to a trio. Each entry should consist of either
#'   one or three numbers. For one number, the locus used for calculating the
#'   summary statistics is locus in the provided data that corresponds to the
#'   number. If three numbers are provided, the locus for calculation is created
#'   by combining the corresponding three loci from the given data.
#' @param ... Additional arguments that will be passed to
#'   \code{\link{as.segsites}}.
#' @export
#' @examples
#' segsites <- create_segsites(matrix(c(1, 0, 0,
#'                                      1, 1, 0,
#'                                      0, 0, 1), 3, 3, TRUE),
#'                             c(.1, .3, .5))
#' model <- coal_model(3, 1) +
#'   sumstat_sfs() +
#'   sumstat_nucleotide_div() +
#'   sumstat_mcmf()
#' sumstats <- calc_sumstats_from_data(model, list(segsites))
#' print(sumstats)
calc_sumstats_from_data <- function(model,
                                    segsites_list = NULL,
                                    tree_list = NULL,
                                    trios = NULL,
                                    ...) {

  assert_that(is.model(model))

  if (!is.null(trios)) {
    assert_that(is.list(trios))
    assert_that(!is.null(segsites_list))
    if (!is.null(tree_list)) {
      stop("Using tree_list with trios is currently not supported")
    }

    segsites_list <- lapply(trios, function(trio) {
      assert_that(is.numeric(trio))
      if (length(trio) == 1) return(segsites_list[[trio[1]]])
      assert_that(length(trio) == 3)
      create_locus_trio(list(segsites_list[[trio[1]]]),
                        list(segsites_list[[trio[2]]]),
                        list(segsites_list[[trio[3]]]))[[1]]
    })
  }

  if (!is.null(segsites_list)) {
    segsites_list <- as.segsites(segsites_list, ...)
    if (!all(vapply(segsites_list, is_segsites, logical(1)))) {
      stop("Incorrect or missing data in list of segregating sites")
    }
    assert_that(length(segsites_list) == get_locus_number(model))
  }

  if (!is.null(tree_list)) {
    stop("Trees are currently not supported")
    #assert_that(is.list(tree_list))
    #assert_that(length(tree_list) == get_locus_number(model))
  }

  calc_sumstats(model, segsites_list, tree_list)
}


requires_segsites <- function(model) {
  any(sapply(get_summary_statistics(model), function(x) x$requires_segsites()))
}


requires_trees <- function(model) {
  any(sapply(get_summary_statistics(model), function(x) x$requires_trees()))
}


requires_files <- function(model) {
  any(sapply(get_summary_statistics(model), function(x) x$requires_files()))
}
