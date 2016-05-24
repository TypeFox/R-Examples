listTypes <- function (method, class = "default", strict = FALSE)
{
	## List all custom functions for a method and for a given class
	## For instance, a custom view is a function as 'view_<customview>.class'

	makeList <- function (method, class) {
		Pat <- paste("^", method, "_([^.]+)\\.", class, "$", sep = "")
		return(sub(Pat, "\\1",
			apropos(Pat, ignore.case = FALSE, mode = "function")))
	}

	method <- as.character(method[1])
	class <- as.character(class)
	List <- makeList(method, class[1])
	if (!isTRUE(strict)) {
		## Also include views for inherited classes
		if (l <- length(class) > 1)
			for (i in 2:l)
				List <- c(List, makeList(method, class[i]))
		## Also include default views
		if (class != "default")
			List <- c(List, makeList(method, "default"))
		## Sort the final list and eliminate duplicated items
		List <- sort(unique(List))
	}
	List
}
