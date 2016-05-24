#' Creates all group sizes by reducing one group in all rows of grpsizes.
#'
#' Used for generating all group size combinations for one specific total size
#' iteratively, starting from grpsizes with one row containing original group
#' sizes.
#'
#' @param grpsizes  A data.table with the columns containing the group names,
#'                  and the rows containing a particular setup of group sizes.
#'                  All rows are expected to have the same sum (not checked).
#'
#' @param grpnames  The group names (specified because the table can have other
#'                  columns as well).
#'
#' @return A data.table with the same format as grpsizes, containing all
#' possible group setups totaling to one less than the total in grpsizes.
decrease_group_sizes <- function(grpsizes, grpnames) {
    d <- rbindlist(lapply(grpnames, function(col) {
        dg <- data.table::copy(grpsizes)
        dg[, col := get(col) - 1, with = FALSE]
        dg <- dg[eval(as.name(col)) > 0]
    }))
    data.table::setkeyv(d, grpnames)
    unique(d)
}


#' Creates Cartesian product of iterators.
#'
#' @param initializers  A list of initializer functions (with no arguments)
#'                      for iterators.
#'
#' @param get_next      A function for retrieving next item for an iterator
#'                      argument; it assumes that the iterator returns NULL
#'                      when finished.
#'
#' @param sspace        elements to be used (a list of vectors)
#'
#' @return A function that returns list of values, and stops with
#'   "StopIteration" message when finished, so that it can be used with the
#'   iterators::iter() function to create an iterator that works with foreach.
create_Cartesian_iterable <- function(initializers, get_next, sspace) {
    stopifnot(length(initializers) == length(sspace))
    values <- NULL
    iterators <- sapply(initializers, function(fn) fn(), simplify = FALSE)
    function() {
        if (is.null(values)) {
            values <<- lapply(seq_along(iterators), function(pos)
                sspace[[pos]][get_next(iterators[[pos]])])
        } else {
            pos <- length(iterators)
            repeat {
                index <- get_next(iterators[[pos]])
                if (!is.null(index)) {
                    values[[pos]] <<- sspace[[pos]][index]
                    inc(pos) <- 1
                    if (pos > length(iterators))
                        break
                } else {
                    iterators[[pos]] <<- initializers[[pos]]()
                    dec(pos) <- 1
                    if (!pos) {
                        values <<- NULL
                        stop("StopIteration")
                    }
                }
            }
        }
        values
    }
}


#' Combines matched group candidate with current best one.
#'
#' @param candidate  A list containing ratio (a number) and ind (a vector).
#'
#' @inheritParams check_subspaces_for_group_size_setup
#'
#' @return A structure like best containing the ind vectors with the highest
#' ratio.
combine_ratios <- function(best, candidate) {
    if (candidate$ratio && candidate$ratio >= best$ratio) {
        if (candidate$ratio > best$ratio) {
            best$ratio <- candidate$ratio
            best$inds <- list()
        }
        best$inds <- c(best$inds, list(candidate$ind))
    }
    best
}


#' Searches over all possible subspaces for specified group size setup.
#'
#' @param best       The best matched groups so far together with its
#'                   p-value / thresh ratio; a list containing ratio and inds
#'                   (a list of subject index vectors).
#'
#' @param grpsize_setup  A set of group sizes as a data.table row (also a list).
#'
#' @inheritParams ldamatch
#'
#' @inheritParams search_exhaustive
#'
#' @return The best
#'
#' @importFrom iterators iter
check_subspaces_for_group_size_setup <- function(
        best, grpsize_setup, sspace, condition, covariates, halting_test, thresh,
        print_info) {
    ci <- create_Cartesian_iterable(
        sapply(names(sspace), function(cond) {
            function() iterpc::iterpc(
                length(sspace[[cond]]), grpsize_setup[[cond]])
        }, simplify = FALSE), iterpc::getnext, sspace)
    if (print_info) {
        Cartesian_size <- prod(vapply(
            get("iterators", environment(ci)), iterpc::getlength, 0))
        cat("Size of Cartesian product:", Cartesian_size, "\n")
        start_time <- proc.time()
    }
    # calculate halting test values for batch
    values <- NULL  # just to eliminate the "NOTE" given by R CMD check
    best <- foreach::foreach(
            values = iterators::iter(ci),
            .combine = combine_ratios, .init = best, .inorder = FALSE) %dopar% {
        ind <- unlist(values)
        ratio <- halting_test(condition[ind], covariates[ind, , drop = FALSE], thresh)
        list(ratio = ratio, ind = ind)
    }
    if (print_info) {
        elapsed_time <- (proc.time() - start_time)[['elapsed']]
        cat("Number of cases processed per second:",
            Cartesian_size / elapsed_time, "\n")
    }
    best
}


#' Searches the space backwards, prefering more subjects and better ratios.
#' While the search is done in parallel, the search space is enormous and so
#' it can be very slow in the worst case. It is perhaps most useful as a tool
#' to study other matching procedures.
#'
#' @param sspace  An ordered subject subspace: a list of vectors,
#' with one vector per group containing the corresponding subject indices.
#'
#' @inheritParams ldamatch
#'
#' @return A list of logical vector(s) for best set(s) of subjects to be kept.
#'
#' @import data.table
#' @import foreach
#' @importFrom entropy KL.plugin
#' @importFrom iterpc iterpc
search_exhaustive <- function(sspace, condition, covariates,
                              halting_test, thresh, props,
                              print_info) {
    # Finds best p-value / threshold ratio and the corresponding subsets:
    # iterates over all group sizes, by decreasing groups with original size.
    best <- list(ratio = 0, inds = list())
    grpsizes <- data.table::data.table(t(vapply(sspace, length, 0)))
    grpnames <- names(sspace)
    while (!best$ratio) {
        grpsizes <- decrease_group_sizes(grpsizes, grpnames)
        if (nrow(grpsizes) == 0)
            break
        if (print_info)
            cat("Created", nrow(grpsizes), "group size configurations",
                "each with a total size of",
                sum(grpsizes[1, grpnames, with = FALSE]), "\n")
        # Orders rows by similarity to ratio in props.
        grpsizes[
            , KL_diverg := vapply(seq_len(nrow(.SD)), function(row)
                entropy::KL.plugin(prop.table(.SD[row, ]), props), 0.0),
            .SDcols = grpnames]
        data.table::setorder(grpsizes, KL_diverg)
        # Goes over rows of table with group sizes and
        # finds best suitable subsets, comparing those for group sizes with the
        # same KL divergence only.
        grpsizes_row <- 1
        repeat {
            KL_diverg <- grpsizes[grpsizes_row, KL_diverg]
            repeat {
                if (print_info)
                    cat(paste(names(grpsizes), grpsizes[grpsizes_row], sep = ": "), "\n")
                best <- check_subspaces_for_group_size_setup(
                    best, grpsizes[grpsizes_row, ], sspace, condition, covariates,
                    halting_test, thresh, print_info)
                inc(grpsizes_row) <- 1
                if (grpsizes_row > nrow(grpsizes) ||
                        grpsizes[grpsizes_row, KL_diverg] != KL_diverg)
                    break
            }
            if (best$ratio || grpsizes_row > nrow(grpsizes))
                break
        }
    }
    if (!best$ratio)
        stop("No subspace found for specified threshold")
    lapply(best$inds, function(ind) seq_along(condition) %in% ind)
}
