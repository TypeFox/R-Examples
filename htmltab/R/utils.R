is_empty <- function(x) length(x) == 0

equal_length <- function(x) {
  x.length <- sapply(x, length)
  length(unique(x.length)) == 1
}

rm_str_white <- function(el) gsub("^\\s+|\\s+$", "", el)

check_correct <- function(nodeset){
  names <- sapply(1:length(nodeset), function(index) XML::xmlName(nodeset[[index]]))
  if(!all(names == "tr")) stop("You must pass header/body information that identifies row elements (tr)", call. = FALSE)
}

#' Wrapper around if stop logic
#' @noRd
ifstop <- function(cond, mess, ...){

  cond <- eval(quote(cond))

  if(isTRUE(cond)){
    stop(mess, call. = F)
  }
}

#' Is str a URL?
#' @noRd
is_url <- function(str){
  !grepl("<?>", str)
}


#' Assert a specific tag in an XML node
#'
#' @param table.Node the table node
#' @param tag a character string for the tag name to be matched
#' @return logical value indicating whether tag is present in table code
#' @noRd

has_tag <- function(table.Node, tag) {
  x <- XML::xpathSApply(table.Node, tag)
  if(length(x) > 0){TRUE} else{FALSE}
}
