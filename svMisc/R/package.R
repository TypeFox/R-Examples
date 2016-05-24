package <- function (..., warn = TRUE)
{
	## A multiple require proceeding as silently as possible
	## Suppress packages messages as much as possible
	owarn <- getOption("warn")
	options(warn = -1)  # Suppress warnings
	on.exit(options(warn = owarn))
	args <- unlist(list(...))
	l <- length(args)
	check <- rep(TRUE, l)
	if (l > 0)
		for (i in 1:l)
			check[i] <- suppressPackageStartupMessages(require(args[i],
				quietly = TRUE, character.only = TRUE, warn.conflicts = FALSE))
	if (!all(check) && isTRUE(warn)) {
		bads <- args[!check]
		options(warn = owarn)
		if (length(bads) == 1) {
			warning("Unable to load package ", bads, "!\n")
		} else {
			warning("Unable to load package(s): ",
				paste(bads, collapse = ", "), "!\n")
		}
	}
	invisible(check)
}
