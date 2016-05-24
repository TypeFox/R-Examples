#' Split a String Into Run Chunks
#' 
#' Splits a string into a vector of runs.
#' 
#' @param x A string.
#' @return Returns a list of vectors. 
#' @seealso \code{\link[qdapTools]{loc_split}},
#' \code{\link[qdapTools]{split_vector}}
#' @author Robert Reed and Tyler Rinker <tyler.rinker@@gmail.com>.
#' @references \url{http://stackoverflow.com/a/24319217/1000343} 
#' @export
#' @examples
#' run_split(c("122333444455555666666", NA, "abbcccddddeeeeeffffff"))
run_split <- function(x) {
    strsplit(x, "(?<=(\\w))(?!\\1)", perl = TRUE)
}