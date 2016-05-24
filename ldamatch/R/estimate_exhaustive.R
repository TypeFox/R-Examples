#' Returns human readable format for number of seconds.
#'
#' @param seconds       The number of seconds to convert to human-readable form.
#'
#' @param num_decimals  The number of decimals to print in the output.
#'
#' @return A string containing "<number> seconds/minutes/hours/days/years".
get_human_readable <- function(seconds, num_decimals = 3) {
    value <- seconds
    measurement <- "second"
    for (convert in list(list(60, "minute"),
                        list(60, "hour"),
                        list(24, "day"),
                        list(7, "week"),
                        list(365/(12*7), "month"),
                        list(12, "year"))) {
        if (round(value, num_decimals) >= convert[[1]]) {
            value <- value / convert[[1]]
            measurement <- convert[[2]]
        } else {
            break
        }
    }
    value <- round(value, num_decimals)
    paste(value, if (value == 1) measurement else paste0(measurement, "s"))
}


#' Estimates the maximum number of cases to be checked during exhaustive search.
#'
#' @param min_preserved     Assumes that at least a total of this many subjects
#'                          will be preserved.
#'
#' @param cases_per_second  Assumes that this number of cases are checked out
#'                          per second, for estimating the time it takes to run
#'                          the exhaustive search; default: 100.
#'
#' @param print_info        If TRUE, prints partial calculations as well for
#'                          the number of cases and estimated time when removing
#'                          1, 2, ... subjects.
#'
#' @inheritParams ldamatch
#'
#' @return The maximum number of cases.
#'
#' @examples
#' estimate_exhaustive(58, as.factor(c(rep('ALN', 25), rep('TD', 44))))
#' estimate_exhaustive(84, as.factor(c(rep('ASD', 51), rep('TD', 44))))
#'
#' @import data.table
#' @importFrom iterpc iterpc
#'
#' @export
estimate_exhaustive <- function(
        min_preserved, condition, cases_per_second = 100,
        print_info = get("PRINT_INFO", .ldamatch_globals)) {
    stopifnot(is.factor(condition),
              min_preserved >= length(levels(condition)))
    sspace <- split(seq_along(condition), condition)
    grpnames <- names(sspace)
    num_cases <- 0
    grpsizes <- data.table::data.table(t(vapply(sspace, length, 0)))
    while (min_preserved < sum(grpsizes[1, names(sspace), with = FALSE])) {
        grpsizes <- decrease_group_sizes(grpsizes, grpnames)
        for (grpsizes_row in 1:nrow(grpsizes)) {
            inc(num_cases) <- prod(sapply(names(sspace), function(cond)
                iterpc::getlength(iterpc::iterpc(
                    length(sspace[[cond]]), grpsizes[[grpsizes_row, cond]]))))
        }
        if (print_info) {
            cat("If", sum(grpsizes[1, names(sspace), with = FALSE]), "of",
                length(condition), "kept: at most", num_cases, "cases. ")
            cat("If ", cases_per_second, " cases per second evaluated: ",
                get_human_readable(num_cases / cases_per_second), ".\n", sep = "")
        }
    }
    num_cases
}
