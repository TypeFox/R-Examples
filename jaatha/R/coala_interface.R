#' Use a coala model in Jaatha
#' 
#' This creates a Jaatha model from a coala model. Simulation for this model
#' model are conducted via the \code{simulate} function for the coala model.
#' The parameters that are
#' estimated must be specified via \code{\link[coala]{par_range}} and the
#' model must not have any other named parameters. Summary statistics present 
#' in the coala model are used in Jaatha.
#' 
#' @param x The coala model
#' @param jsfs_summary The way the Joint Site Frquency Spectrum (JSFS) 
#'   is further summarized. Can be \code{sums} (default), \code{none} or 
#'   \code{"smoothing"}. For \code{sums}, 23 different areas of the JSFS
#'   are summed up, and the sums are used as indepented Poission statistcs, 
#'   for \code{none}, all entries are used as indepented Possion statistics.
#'   The value \code{smooth} is experimental so far and should not be used.
#'   This option has no effect if the JSFS used as summary statistic in the
#'   coala model.
#' @param four_gamete_breaks Quantiles of the real data that will be used as 
#'   breaks for binning the Four Gamete test based statistic if present in the 
#'   model.
#' @param mcmf_breaks Quantiles of the real data that will be used as breaks
#'   for binning the MCMF statistic if present in the model.
#' @param jsfs_part Partitions used for the summarizing the JSFS. This is only 
#'   used if \code{jsfs_summary} is "sums". Is used as the \code{part} argument
#'   of \code{\link{coarsen_jsfs}}. Please go there for an explanation.
#' @param jsfs_part_hi Same as \code{jsfs_part}, but used as \code{part_hi} 
#'   argument in \code{\link{coarsen_jsfs}}.
#' @inheritParams create_jaatha_model
#' @export
#' @export create_jaatha_model.coalmodel
create_jaatha_model.coalmodel <- function(x, 
                                          jsfs_summary = c("sums",
                                                           "none",
                                                           "smooth"),
                                          four_gamete_breaks = c(.2, .5),
                                          mcmf_breaks = c(.5, .7, .9),
                                          jsfs_part = c(1, 3),
                                          jsfs_part_hi = c(1, 3),
                                          ...,
                                          scaling_factor = 1,
                                          test = TRUE) {
  
  require_package("coala")

  if (length(jsfs_summary) > 1) jsfs_summary <- jsfs_summary[1]
  
  sim_func <- function(pars) stats::simulate(x, pars = pars)
  
  # create parameter ranges
  par_table <- coala::get_parameter_table(x)
  par_ranges <- as.matrix(par_table[, -1])
  rownames(par_ranges) <- par_table$name

  # create summary statisics
  sum_stats <- convert_coala_sumstats(x, jsfs_summary,
                                      four_gamete_breaks, mcmf_breaks,
                                      jsfs_part, jsfs_part_hi)
  
  create_jaatha_model.function(sim_func, par_ranges, sum_stats, 
                               test = test)
}


convert_coala_sumstats <- function(coala_model, jsfs_summary = "sums",
                                   four_gamete_breaks, mcmf_breaks,
                                   jsfs_part, jsfs_part_hi) {
  
  require_package("coala")
  assert_that(is.string(jsfs_summary))
  
  lapply(coala::get_summary_statistics(coala_model), function(stat) {
    name <- stat$get_name()
    
    # --- JSFS Summary Statistic ------------------------------------
    if (inherits(stat, "stat_jsfs")) {
      if (jsfs_summary == "sums") {
        return(create_jaatha_stat(name, function(x, opts) {
          coarsen_jsfs(x[[name]], jsfs_part, jsfs_part_hi)
        }))
      } else if (jsfs_summary == "none") {
        return(create_jaatha_stat(name, function(x, opts) {
          as.vector(x[[name]])[-c(1, prod(dim(x[[name]])))]
        }))
      } else if (jsfs_summary == "smooth") {
        stop("Smoothing is not suppored right now")
      }
    }
    
    # --- SFS Summary Statistic ------------------------------------
    if (inherits(stat, "stat_sfs")) {
      return(create_jaatha_stat(name, function(x, opts) x[[name]]))
    }
    
    # --- Four Gamete Summary Statistic -----------------------------
    if (inherits(stat, "stat_four_gamete")) {
      return(create_jaatha_stat(name, function(x, opts) {
        x[[name]][, c(1, 2, 6), drop = FALSE]
      }, poisson = FALSE, breaks = four_gamete_breaks))
    }
    
    # --- OmegaPrime Summary Statistic ----------------------------------
    if (inherits(stat, "stat_omega_prime") || inherits(stat, "stat_mcmf")) {
      return(create_jaatha_stat(name, function(x, opts) x[[name]],
                                poisson = FALSE, 
                                breaks = mcmf_breaks))
    }
    
    warning("Summary statistic '", name, "' is not supported. Ignoring it.")
    NULL
  })
}


multi_index_range <- function(d, p) {
  ## d are dimensions of an array A, and p is a matrix of numbers. Then this
  ## function returns a vector v such that
  ## A[p[1,1]:p[1:2],p[2,1]:p[2:2],...] consists of the same values as A[v],
  ## even though no necessarily in the same order.
  N <- nrow(p)
  v <- p[N, 1]:p[N, 2]
  if (N > 1) {
    for (n in (N - 1):1) {
      v <- as.vector(outer( (v - 1) * d[n], p[n, 1]:p[n, 2], "+"))
    }
  }
  
  v
}


#' Divides the joint site frequency spectrum (jsfs) into blocks
#' and returns the sum of the jsfs entries for each block.
#' 
#' ja is the jsfs, part a list of vectors specifying for each dimension
#' how ja should be partitioned. If part_hi!=NULL, it is a list spefifying
#' how ja is to be paritioned on the higher end of each dimension.  if
#' part or part_hi is not a list, it is turned into a list of the same
#' length as dim(ja), in which each entry is the original part or part_hi
#' e.g. 2,7,9 partitions into 1:2, 3:7, 8:9, 9:N For example, with
#' part=c(1,3) and part_hi=c(1,3) we get the classical jaatha summary
#' statistics. Note, however, that the order in which they appear will be
#' different than in the original jaatha package.
#' 
#' @param ja an array containing the joint site frequency spectrum
#' @param part a vector of integers or a list of vectors of integers. If
#'   it is a list, the vector part[[i]] specifies that the \eqn{i}-th dimension
#'   of \code{ja} should be partitioned into \code{1:(part[[i]][1]-1)},
#'   \code{part[[i]][1]:(part[[i]][2]-1)}, and so on. If \code{part} is a
#'   vector, it will be used for all dimensions.
#' @param part_hi NULL or a vector of integers or a list of vector of integers
#'    indicating the partioning at the higher end of each dimension. This means,
#'    if it is a list, the values in the vector \code{dim(ja)[i]-part_hi[[i]]}
#'    will be appended to the end of \code{part[[i]]}. If \code{part_hi} is a
#'    single vector, it will be used for all dimensions. Thus, with the
#'    combination of part=c(1,3) and part_hi=c(1,3), the classical jaatha summary
#'    statistics, plus the two values \code{ja[0]} and
#'    \cite{ja[length(ja)]}. Note that the order in which they appear will
#'    however be different than in the original jaatha summary statistics.
#' @return vector of numbers, which are the sums over the blocks of the jsfs
#'    for all combinations of partitions
#' @author Dirk Metzler & Paul Staab
#' @references A. Tellier, P. Pfaffelhuber, B. Haubold, L. Naduvilezhath,
#'   L. E. Rose, T. Staedler, W. Stephan, and D. Metzler (2011) Estimating
#'   parameters of speciation models based on refined summaries of the joint
#'   site-frequency spectrum. PLoS One 6(5): e18155
coarsen_jsfs <- function(ja, part, part_hi = NULL) {
  d <- dim(ja)
  n <- length(d)
  if (!is.list(part)) part <- rep(list(part), n)
  if (!is.null(part_hi)) {
    if (!is.list(part_hi)) part_hi <- rep(list(part_hi), n)
    for (i in 1:n) {
      upper <- sort(d[i] - part_hi[[i]])
      if (utils::tail(part[[i]], 1) >= upper[1]) {
        stop(paste("part and part_hi incompatible in dim", i))
      }
      part[[i]] <- c(part[[i]], upper)
    }
  }
  
  for (i in 1:n) {
    part[[i]] <- c(0, part[[i]], dim(ja)[i])
  }
  
  z <- numeric(length = prod(vapply(part, length, numeric(1)) - 1))
  combinations <- 
    expand.grid(lapply(vapply(part, length, numeric(1)) - 1, ":", 1))[length(z):1, ] #nolint
  
  for (i in 1:length(z)) {
    comb <- combinations[i, ]
    p <- matrix(NA, ncol = 2, nrow = n)
    for (j in 1:n) {
      p[j, 1] <- part[[j]][comb[[j]]] + 1
      p[j, 2] <- part[[j]][comb[[j]] + 1]
    }
    z[i] <- sum(ja[multi_index_range(d, p)])
  }
  
  if (all(vapply(part, function(x) any(x == 1), logical(1)))) z <- z[-1]
  if (all(mapply(function(x, y) any(x == y), part, d - 1))) z <- z[-length(z)]
  
  z
}
