

#' \Sexpr{tools:::Rd_package_title("TSTr")}
#' 
#' \Sexpr{tools:::Rd_package_description("TSTr")}
#' 
#' 
#' This package can be used to create a ternary search tree, more space efficient compared to standard prefix trees. Common applications for ternary search trees include spell-checking and auto-completion.
#' 
#' @name TSTr-package
#' @aliases TSTr-package TSTr
#' @docType package
#' @author
#' \Sexpr{tools:::Rd_package_author("TSTr")}
#' 
#' Maintainer:
#' \Sexpr{tools:::Rd_package_maintainer("TSTr")}
#' @seealso \code{\link{newTree}}
#' @references \url{https://en.wikipedia.org/wiki/Ternary_search_tree}
#' @keywords package
#' 
NULL





#' 10001 most frequent English words
#' 
#' A character vector containing 10001 frequency-ordered English words. All words have 3 or more letters.
#' 
#' @name XMIwords
#' @docType data
#' @format The format is: chr [1:10001] "the" "and" "that" "for" "you" "with"
#' "was" "this" "have" "but" "are" "not" "from" ...
#' @references Extracted from the English HC Corpora.
#' @keywords datasets
#' @examples
#' 
#' data(XMIwords)
#' str(XMIwords)
#' 
NULL



