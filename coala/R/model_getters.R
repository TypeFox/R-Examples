#' Getters for coalescent models
#'
#' @param model The coalescent model from which aspects are returned
#'
#' @export
#' @keywords internal
#' @author Paul Staab
get_features <- function(model) model$features


#' @export
#' @describeIn get_features Returns the ranged parameters of a model as a
#'   data.frame
get_parameter_table <- function(model) {
  stopifnot(is.model(model))

  par_table <- read_cache(model, "par_table")
  if (is.null(par_table)) {
    if (!all(sapply(get_parameter(model), is.ranged_par))) {
      stop("Can not create a parameter table with non-ranged pars in model")
    }
    if (length(get_parameter(model)) == 0) {
      par_table <- (data.frame(name = character(),
                               lower.range = numeric(),
                               upper.range = numeric(),
                               stringsAsFactors = FALSE))
    } else {
      par_table <- do.call(rbind, lapply(get_parameter(model), function(par) {
        data.frame(name = par$get_name(),
                   lower.range = par$get_range()[1],
                   upper.range = par$get_range()[2],
                   stringsAsFactors = FALSE)
      }))
    }
    cache(model, "par_table", par_table)
  }

  par_table
}


#' @export
#' @describeIn get_features Returns the ranged parameters of a model
get_parameter <- function(model) {
  stopifnot(is.model(model))
  model$parameter
}


#' @param locus The number of the locus.
#' @param total If \code{FALSE}, the length of loci in a trio will be reported
#'   individually. If \code{TRUE} the sum of the loci"s length will be reported.
#'   This does not affect non-trio loci.
#' @param group The group of loci.
#'
#' @describeIn get_features Returns the length of the loci in a locus group
#' @export
get_locus_length <- function(model, locus = NULL, group = NULL, total = TRUE) {
  llm <- get_locus_length_matrix(model)
  if (is.null(locus) & is.null(group)) {
    group <- 1:nrow(llm)
  }

  # Group and locus are identical for ilv models
  if (!is.null(group) && has_variation(model)) {
    locus <- group
  }

  # Determine to which group the locus belongs
  if (!is.null(locus)) {
    group <- get_locus_group(model, locus)
  }

  if (total) return(rowSums(llm[group, 1:5, drop = FALSE]))
  if (sum(llm[group , c(1:2, 4:5)]) == 0) return(llm[group , 3])
  llm[group, 1:5]
}


get_locus_group <- function(model, locus) {
  llm <- get_locus_length_matrix(model)
  min(which(cumsum(llm[ , "number"]) >= locus))
}


get_locus_group_number <- function(model) {
  nrow(get_locus_length_matrix(model))
}


#' @describeIn get_features Returns a vector of populations in the model
#' @export
get_populations <- function(model) {
  feat <- model$features[vapply(model$features, is_feat_sample, logical(1))]
  pop_number <- max(vapply(feat, function(x) length(x$get_sizes()), numeric(1)))
  assert_that(pop_number > 0)
  1:pop_number
}


#' @describeIn get_features Returns a matrix with detailed length
#' information about the loci in the model.
#' @export
get_locus_length_matrix <- function(model) {
  llm <- read_cache(model, "llm")
  if (is.null(llm)) {
    assert_that(length(model$loci) >= 0)
    llm <- do.call(rbind, lapply(model$loci, function(l) {
             number <- ifelse(l$get_number() > 1,
                              round(l$get_number() / model$scaling_factor),
                              l$get_number())
             c(l$get_length(TRUE), number = number)
           }))

    cache(model, "llm", llm)
  }
  llm
}


#' @describeIn get_features Returns the number of loci in a locus group
#' @param ignore_variation For internal use. Will likely be removed soon.
#' @export
get_locus_number <- function(model, group = NA, ignore_variation = FALSE) {
  numbers <- get_locus_length_matrix(model)[ , "number"]
  if (is.na(group)) return(sum(numbers))
  if (has_variation(model) && !ignore_variation) return(1)
  numbers[group]
}


#' @describeIn get_features Returns the index of the individuals of one
#'   population. Ignores outgroups, so that it can be used for indexing
#'   segregating sites.
#' @param zero_indexed If true, the names of the populations are started from
#'   0 instead of from 1.
#' @param haploids If \code{TRUE}, the function always returns all haploids
#'   from the population, even if the model is polyploid.
#' @export
get_population_individuals <- function(model, pop,
                                      zero_indexed = FALSE,
                                      haploids = TRUE) {

  sample_size <- get_sample_size(model)
  outgroup <- get_outgroup(model)

  if (!is.na(outgroup)) {
    if (pop == outgroup) stop("Calculating summary statistics for the outgroup")
    sample_size[outgroup] <- 0
  }

  if (!haploids) sample_size <- sample_size / get_ploidy(model)
  if (pop == "all") return(1:sum(sample_size))

  if (!pop %in% get_populations(model)) stop("Invalid population")
  from <- cumsum(c(0, sample_size)) + 1
  to <- cumsum(sample_size)
  from[pop]:to[pop]
}


get_par_names <- function(model, without_priors=FALSE) {
  param <- get_parameter(model)
  if (length(param) == 0) return(character(0))
  if (without_priors) {
    param <- param[!vapply(param, is.prior_par, numeric(1))]
  }
  sapply(param, function(par) par$get_name())
}

get_loci <- function(model) model$loci

get_cmd <- function(model) {
  suppressWarnings(prog <- select_simprog(model))
  if (is.null(prog)) stop("No suitable simulator found.")
  prog$get_cmd(model)
}
