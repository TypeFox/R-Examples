#'
#' A data frame of ballpark IDs
#'
#' This function returns a two-column data frame of ballpark IDs
#' along with current stadium name
#'
#' @export
#'
#' @examples getParkIDs()
#'

getParkIDs <- function() {
    u <- "http://www.retrosheet.org/parkcode.txt"
    nm <- scan(u, nlines = 1L, nmax = 2L, what = "", sep = ",", quiet = TRUE)
    what <- setNames(rep_len(list(""), length(nm)), nm)
    scn <- scan(u, skip = 1L, what = what, sep = ",", flush = TRUE,
        quote = ",", quiet = TRUE)
    attr(scn, "row.names") <- .set_row_names(length(scn[[1L]]))
    class(scn) <- "data.frame"
    scn
}
