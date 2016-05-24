##' Extract pollen or other proxy counts from data objects and returns them in a useful format.
##'
##' Methods are available for "download" and "download_list" objects.
##'
##' @title Access proxy count data
##'
##' @param obj an R object from which counts are to be extracted.
##' @param ... arugments passed to other methods.
##' @return Either a data frame of counts or a list of such objects.
##'
##' @author Gavin Simpson
##'
##' @export
##' @rdname counts
##'
##' @examples
##' \dontrun{
##' marion <- get_site('Marion Lake%')
##' louise <- get_site('Louise Pond%')
##' western.sites <- rbind(marion, louise)
##' western.data  <- get_dataset(western.sites)
##'
##' western.dl <- get_download(western.data)
##' western.cnt <- counts(western.dl)
##' sapply(western.cnt, dim)
##' marion.cnt<- counts(western.dl[[1]])
##' dim(marion.cnt)
##' }
`counts` <- function(obj, ...) {
    UseMethod("counts")
}

##' @export
##' @rdname counts
`counts.download` <- function(obj, ...) {
    ret <- as.data.frame(obj$counts)
    class(ret) <- c("neo_counts", "data.frame")
    ret
}

##' @export
##' @rdname counts
`counts.download_list` <- function(obj, ...) {
    ret <- lapply(obj, '[[', 'counts')
    ret <- lapply(ret, as.data.frame)
    class(ret) <- c("neo_counts_list", "list")
    ret
}
