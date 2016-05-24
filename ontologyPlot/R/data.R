#' Capitalise words in character vector
#'
#' @param x Character vector 
#' @return Character vector
#' @examples
#' simple_cap(c("a simple test", "Another-test"))
#' @export
simple_cap <- function(x) {
	s <- strsplit(x, " ")
	sapply(
		s,
		function(x) paste(
			toupper(substring(x, 1,1)), 
			substring(x, 2),
			sep="", collapse=" "
		)
	)
}

#' Get human readable, shortened (where possible) ontological term names
#'
#' @template ontology
#' @template terms
#' @return Character vector
#' @examples
#' library(ontologyIndex)
#' data(hpo)
#' get_shortened_names(hpo, c("HP:0001873", "HP:0011877"))
#' @export
#' @import magrittr
get_shortened_names <- function(ontology, terms) gsub(
	"Impaired |(Abnormality of (the )?)|(Abnormal )", 
	"", 
	ontology$name[terms]
) %>% 
(function (x) gsub("^\\s+|\\s+$", "", x)) %>%
sapply(simple_cap)
