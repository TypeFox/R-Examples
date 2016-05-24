###############################################################################
#
# webvis: An R package to create web visualizations.
# author: Shane Conway <shane.conway@gmail.com>
#
# This is released under a BSD license.
#
# Documentation was created using roxygen:
# roxygenize('webvis', roxygen.dir='webvis', copy.package=FALSE, unlink.target=FALSE)
#
###############################################################################

#' Checks whether something is a webvis object.
#'
#' \code{is.webvis} Checks whether something is a webvis object.
#'
#' @param x The object to be checked.
#' @return A logical value.
#' @author Shane Conway \email{shane.conway@@gmail.com}
is.webvis <- function(x) class(x) == "webvis"

#' Checks whether something is a webvis.flat object.
#'
#' \code{is.webvis.flat} Checks whether something is a webvis.flat object.
#'
#' @param x The object to be checked.
#' @return A logical value.
#' @author Shane Conway \email{shane.conway@@gmail.com}
is.webvis.flat <- function(x) class(x) == "webvis.flat"

#' Checks whether something is a webvis.param object.
#'
#' \code{is.webvis.param} Checks whether something is a webvis.param object.
#'
#' @param x The object to be checked.
#' @return A logical value.
#' @author Shane Conway \email{shane.conway@@gmail.com}
is.webvis.param <- function(x) class(x) == "webvis.param"

#' Checks whether something exists.
#'
#' \code{esse} Checks whether something exists.
#'
#' @param x The object to be checked.
#' @return A logical value.
#' @author Shane Conway \email{shane.conway@@gmail.com}
esse <- function(x) {
	if(missing(x)) 
		return(FALSE)
	#test <- try(length(x))
	#if(class(test) == "try-error") if(all(grepl("is missing", test))) x <- NULL
	#if(!exists(x)) 
	#	return(FALSE)
	if(!length(x)) 
		return(FALSE)
	if(all(is.null(x)))
		return(FALSE)
	if(all(is.na(x)))
		return(FALSE)
	if(is.character(x) && all(x==""))
		return(FALSE)
	return(TRUE)
}

#' A simplied version of paste.
#'
#' \code{collapse} A simplied version of paste. 
#'
#' @param ... Comma separated list of values.
#' @return A string of the values combined.
#' @author Shane Conway \email{shane.conway@@gmail.com}
#' @examples
#' collapse("a", "b", "c")
collapse <- function(...) paste(c(...), collapse="")

#' Return the tail html.
#'
#' \code{getTail} Return the tail html.
#'
#' @param t The tail html.
#' @return HTML.
#' @author Shane Conway \email{shane.conway@@gmail.com}
getTail <- function(t) {
	if(!esse(t)) t <- "</body></html>"
	return(t)
}

#' Return the head html.
#'
#' \code{getHead} Return the head html.
#'
#' @param title The title of the page.
#' @param protovis.path The path (url) for protovis.
#' @return HTML.
#' @author Shane Conway \email{shane.conway@@gmail.com}
getHead <- function(title="", protovis.path=PROTOVIS.PATH) {
	h <- paste("<html>
					<head>
					<title>", title, "</title>
					<link type='text/css' rel='stylesheet' href='ex.css?3.1'/>
					<script type='text/javascript' src='", protovis.path, "'></script>
					</head>
					<body>")
	return(h)
}

#' Checks whether a field exists in a dataset.
#'
#' \code{field.exists} Checks whether a field exists in a dataset.
#'
#' @param field The field name.
#' @param data The data to be checked.
#' @return A logical value.
#' @author Shane Conway \email{shane.conway@@gmail.com}
field.exists <- function(field, data) {
	if(!esse(data)) FALSE
	if(!missing(data)) if(is.data.frame(data)) { if(any(field %in% colnames(data))) return(TRUE) } else { stop("data should be a data.frame")}
	FALSE
}
