#'
#' Files currently available for download
#'
#' A convenience function, returning the base file names of the
#' available downloads for the \code{year} and \code{type} arguments
#' in \code{getRetrosheet}.
#'
#' @return A named list of available single-season Retrosheet event and
#' game-log zip files, and schedule text files. These file names are
#' not intended to be passed to \code{getRetrosheet}, but is simply a
#' fast way to determine if the desired data is available.
#'
#' @examples getFileNames()
#'
#' @importFrom RCurl getCurlHandle
#' @importFrom RCurl getURL
#' @importFrom XML htmlParse
#' @importFrom XML xpathSApply
#' @importFrom XML free
#'
#' @export

getFileNames <- function() {
    path <- c(event = "game.htm", gamelog = "gamelogs/index.html",
        schedule = "schedule/index.html")
    full <- sprintf("http://www.retrosheet.org/%s", path)
    curl <- getCurlHandle()
    docs <- lapply(full, function(x) {
        content <- getURL(x, curl = curl)
        htmlParse(content, asText = TRUE)
    })
    o <- function(pat, doc) {
        fnames  <- xpathSApply(doc,
            path = "(//pre//a | //b/a)/@href", fun = basename)
        grep(pat, fnames, value = TRUE)
    }
    part <- sprintf(c("%seve.zip", "gl%s.zip", "%ssked.txt"), "\\d+")
    res <- setNames(Map(o, part, docs), names(path))
    lapply(docs, free)
    res
}
