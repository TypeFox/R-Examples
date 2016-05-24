#' ldamatch: Multivariate Condition Matching by Backwards Elimination using Fisher's Linear Discriminant
#'
#' Performs group matching by backward elimination using linear discriminant analysis.
#'
#' @docType package
#' @name GUtils
NULL


## Actual matching procedure.

#' Creates ordered subspace of subject candidates for removal, using LDA.
#'
#' @param condition     A factor vector containing condition labels.
#'
#' @param covariates    A vector or columnwise matrix containing
#'                      covariates to match the conditions on.
#'
#' @return An ordered subject subspace: a list of vectors, with one vector per
#' group containing the corresponding subject indices.
#'
#' @importFrom MASS lda
#' @importFrom stats coef
create_subject_subspace_using_LDA <- function(condition, covariates) {
    ## Computes linear projection vector.
    if (ncol(covariates) > 1)
        W <- coef(MASS::lda(condition ~ covariates))[, 1]
    else # Or, just uses the identity projection.
        W <- 1
    projection <- as.vector(covariates %*% W)
    ## Set up search space.
    # Computes means.
    mu <- mean(projection)
    mu.conditions <- tapply(projection, condition, mean)[condition]
    # Locates observations driving condition/covariate correlation(s).
    correlates <- ((projection < mu.conditions) == (mu.conditions < mu))
    # Computes order, filtering out non-correlates.
    ord <- order(projection)
    ord <- ord[correlates[ord]]
    # Splits on condition and reverse for those above.
    sspace <- split(ord, condition[ord])
    for (name in names(sspace))
        if (projection[sspace[[name]][1]] > mu)
            sspace[[name]] <- rev(sspace[[name]])
    sspace
}


#' Creates a matched group via backward selection.
#'
#' @param condition     A factor vector containing condition labels.
#'
#' @param covariates    A vector or columnwise matrix containing
#'                      covariates to match the conditions on.
#'
#' @param halting_test  A function to apply to `covariates` (in matrix form)
#'                      which is TRUE iff the conditions are matched.
#'
#' @param thresh        The statistical threshold to pass onto the
#'                      aforementioned test.
#'
#' @param method        The choice of search method. The "heuristic" method
#'                      deploys the table of desired proportions
#'                      (see below) to structure search. The "montecarlo"
#'                      method randomly generates subspaces of decreasing
#'                      size. The "exhaustive" considers all possible
#'                      subsamples and may be very, very slow.
#'
#' @param props         The desired proportions (percentage) of the sample for
#'                      each condition; if not specified, the (full)
#'                      sample proportions are used. This is used
#'                      for the "heuristic" and "exhaustive" methods.
#'
#' @param replicates    The maximum number of Monte Carlo replications to be
#'                      performed. This is only used for the "montecarlo"
#'                      method.
#'
#' @param print_info    If TRUE, prints summary information on the input and the
#'                      results, as well as progress information for the
#'                      exhaustive search algorithm. Default: TRUE;
#'                      can be changed using set_param("PRINT_INFO", FALSE).
#'
#' @return              A logical vector, TRUE iff row is in the match,
#'                      or a list of such vectors for the exhaustive search.
#'
#' The exhaustive search method uses the foreach package to parallelize
#' computation. To take advantage of this, you must register a cluster.
#' For example, to use all but one of the CPU cores:
#'   doMC::registerDoMC(max(1, parallel::detectCores() - 1))
#' To use sequential processing:
#'   foreach::registerDoSEQ()
#'
#' @export
ldamatch <- function(condition, covariates, halting_test, thresh = .2,
                     method = c("heuristic", "montecarlo", "exhaustive"),
                     props = NULL, replicates = NULL,
                     print_info = get("PRINT_INFO", .ldamatch_globals)) {
    ## Checks arguments and create set their values if missing.
    method <- match.arg(method)
    covariates <- as.matrix(covariates)
    stopifnot(is.factor(condition),
              is.numeric(covariates),
              is.function(halting_test),
              is.numeric(thresh))
    # Checks props argument.
    if (is.null(props)) {
        props <- prop.table(table(condition))
    } else if (method %in% c("heuristic", "exhaustive")) {
        RUnit::checkTrue(sum(props) == 1.0, "sum of props must be 1.0")
        RUnit::checkTrue(length(props) == nlevels(condition),
                         paste("length of props must be", nlevels(condition)))
    } else {
        warning("props parameter ignored for method ", method)
    }
    # Checks replicates argument.
    if (is.null(replicates)) {
        replicates <- get("MC_DEFAULT_REPLICATES", .ldamatch_globals)
    } else if (method %in% c("montecarlo")) {
        RUnit::checkTrue(length(replicates) == 1 && replicates %% 1 == 0,
                         "replicates parameter must be one integer number")
    } else {
        warning("replicates parameter ignored for method ", method)
    }
    # Checks length of container arguments.
    L <- length(condition)
    stopifnot(L == nrow(covariates))

    ## Checks for a "natural match" before setting up search.
    if (halting_test(condition, covariates, thresh))
        return(if (method == "exhaustive") list(rep(TRUE, L)) else rep(TRUE, L))

    ## Search.
    if (print_info) {
        grpsizes <- table(condition)
        cat("Initial group sizes: ",
            paste(names(grpsizes), grpsizes, sep = ": ", collapse = "\t"), "\n")
        cat("Starting", method, "search.\n")
        search_start_time <- proc.time()
    }
    if (method == "heuristic") {
        sspace <- create_subject_subspace_using_LDA(condition, covariates)
        is.in <- search_heuristic(sspace, condition, covariates,
                                  halting_test, thresh, props)
    } else if (method == "montecarlo") {
        sspace <- create_subject_subspace_using_LDA(condition, covariates)
        is.in <- search_montecarlo(sspace, condition, covariates,
                                   halting_test, thresh, replicates)
    } else if (method == "exhaustive") {
        sspace <- split(seq_along(condition), condition)
        is.in <- search_exhaustive(sspace, condition, covariates,
                                   halting_test, thresh, props, print_info)
    }
    if (print_info) {
        total_search_time <- (proc.time() - search_start_time)[['elapsed']]
        cat("Finished", method, "search in", total_search_time, "seconds.\n")
        grpsizes <- table(condition[if (is.list(is.in)) is.in[[1]] else is.in])
        cat("Eventual group sizes:",
            paste(names(grpsizes), grpsizes, sep = ": ", collapse = "\t"), "\n")
        grpremoved <- table(condition) - grpsizes
        cat("Removed subjects:    ",
            paste(names(grpremoved), grpremoved, sep = ": ", collapse = "\t"), "\n")
    }
    is.in
}
