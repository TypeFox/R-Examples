listMethods <- function (f = character(), class = NULL, S3 = TRUE, S4 = TRUE,
mixed = TRUE, filter = getOption("svGUI.methods"))
{
	## Given a function, if it is generic then return a list of its methods
	## or given a class name, return all methods for this class

	## Check argument
	if (!inherits(f, "character"))
		stop("'f' must ba a character string!")

	## List methods for a given class
	if (!is.null(class)) {
		class <- as.character(class)[1]
		res <- list()

		## S3 version
		if (isTRUE(S3)) {
			s3 <- unclass(methods(class = class))
			attr(s3, "info") <- NULL
			## Do we have to filter the methods?
			if (!is.null(filter))
				s3 <- s3[s3 %in% paste(filter, class, sep = ".")]
			res$S3 <- sub(paste(".", class, sep = ""), "", s3)
		}

		## S4 version
		if (isTRUE(S4)) {
			if (is.null(filter)) filter <- character()
			s4 <- capture.output(showMethods(filter, classes = class,
				inherited = FALSE, showEmpty = FALSE))
			## I need to filter this output to get only function names
			res$S4 <- sub("^.*\\: +([^ ]+).*$", "\\1",
				s4[regexpr("Function:", s4) == 1])
		}
		if (isTRUE(mixed)) res <- sort(unique(c(res$S3, res$S4)))
		return(res)

	} else {  # List all methods for one generic function
		## Keep only first item if a vector is provided
		f <- f[1]
		res <- list()

		## S3 version
		if (isTRUE(S3)) {
			## Does the function exists somewhere?
			if (length(findFunction(f, where = .GlobalEnv)) > 0) {
				s3 <-  unclass(suppressWarnings(methods(f)))
				attr(s3, "info") <- NULL
				## Rework this to match presentation for S4 methods
				arg <- names(formals(eval(parse(text =
				paste("getAnywhere(", f, ")", sep = "")))[1]))[1]
				s3 <- sub(paste("^", f, ".", sep = ""), "", s3)
				if (length(s3) == 0 || s3 == "") {
					res$S3 <- character(0)
				} else {
					## Check all possible methods in turn, to verify them
					for (i in 1:length(s3))
					if (inherits(try(getS3method(f, s3[i]), silent = TRUE),
						"try-error")) s3[i] <- ""
					s3 <- s3[s3 != ""]
					if (length(s3) == 0) res$S3 <- character(0) else
						res$S3 <- paste(arg, "=\"", s3, "\"", sep = "")
				}
			} else {  # Not found
				res$S3 <- character(0)
			}
		}

		## S4 version
		if (isTRUE(S4)) {
			## Is it an S4 generic function?
			if (isGeneric(f, where = .GlobalEnv)) {
				s4 <- capture.output(showMethods(f,
					inherited = FALSE, showEmpty = FALSE))
				res$S4 <- s4[-c(1, length(s4))]
			} else {
				res$S4 <- character(0)
			}
		}
		if (isTRUE(mixed)) res <- sort(unique(c(res$S3, res$S4)))
		return(res)
	}
}
