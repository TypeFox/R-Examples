utils::globalVariables(c("bysize", "Mb"))

#' Calculate Sizes of Objects in Workspace
#'
#' Calculate the sizes of all of the objects in the workspace.
#'
#' @param obj Vector of object names. If missing, pull out all object names.
#' @param bysize If \code{TRUE}, sort the objects from smallest to largest.
#'
#' @importFrom plyr arrange
#'
#' @details
#' Calls \code{\link[utils]{object.size}} to get the sizes of a list of objects.
#'
#' @return
#' A data frame with the only column being the size of each object in
#' megabytes (Mb). The row names are the names of the objects.
#'
#' @export
#'
#' @examples
#' print(output <- objsizes())
#' \dontrun{sum(output)}

objsizes <- function(obj, bysize = TRUE) {
    if (missing(obj)) {
      obj <- objects(pos = 1)
    }
    result <- data.frame(name = rep(NA, length(obj)), Mb = rep(0, length(obj)))
    result$name <- obj
    for (i in seq(along = obj)) {
      result[i, 2] <- (utils::object.size(get(obj[i], pos = 1)) / 1024 ^ 2)
    }
    if (bysize == TRUE) {
      result <- plyr::arrange(result, Mb, decreasing = TRUE)
    }
    return(result)
}
