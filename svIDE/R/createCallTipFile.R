createCallTipFile <- function (file = "Rcalltips.txt", pos = 2:length(search()),
field.sep = "=", only.args = FALSE, return.location = FALSE)
{
	## Create a .txt file containing calltips for R functions.
	cat("", file = file) # Create the beginning of the file

	## Get the list of keywords
	keys <- getKeywords(pos = pos)

	## For each keyword, write a line in the file with keyword=calltip
    for (key in keys) {
        ctip <- callTip(key, only.args = only.args)
        if (ctip != "") {
            if (return.location == TRUE) {
				## Get the package from where it is located and append it
				pkg <- sub("^package:", "", find(key, mode = "function"))
				if (length(pkg) > 0 && pkg != ".GlobalEnv")
					pkg <- paste(" [", pkg, "]", sep = "") else pkg <- " []"
            } else pkg <- ""
            cat(key, field.sep, ctip, pkg, "\n", sep = "", file = file,
				append = TRUE)
		}
    }
}
