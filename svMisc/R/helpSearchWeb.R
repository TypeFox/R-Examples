helpSearchWeb <- function (what, type = c("R", "archive", "wiki", "google"),
browse = TRUE, msg = browse, ...)
{
	what <- paste(what, collapse = " ", sep = "")
	what <- gsub(" ", "+", what)
	type <- match.arg(type)
	searchURL <- switch(type,
		"R" = RSiteSearch(what, ...),
		"archive" = paste("http://www.google.com/u/newcastlemaths?q=",
			what, sep = ''),
		"wiki" = paste("http://rwiki.sciviews.org/doku.php?do=search&id=",
			what, sep = ''),
		"google" = paste("http://www.google.com/search?sitesearch=r-project.org&q=",
			what, sep = ''),
		stop("'type' could be only 'R', 'archive', 'wiki' or 'google', currently!"))
	if (type != "R") {
		if (isTRUE(browse)) browseURL(searchURL)
		if (isTRUE(msg)) {
			cat(gettext("A search query has been submitted"), "\n")
			cat(gettext("The results page should open in your browser shortly\n"))
		}
	}
	invisible(searchURL)
}
