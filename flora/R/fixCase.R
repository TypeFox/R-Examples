#' Fix the name case of a taxon
#' 
#' @param x a unit character vector with a taxon
#' @export
#' @examples
#' fixCase("myrcia lingua")
#' fixCase("Myrcia Lingua")
#' fixCase("COFFEA ARABICA")
fixCase <-
function(x) {
    s <- paste(tolower(strsplit(x, " ")[[1]]), collapse = " ")
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}