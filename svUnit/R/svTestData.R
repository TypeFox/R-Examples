is.svTestData <- function (x)
{
	## It this a svTestData object
	return(inherits(x, "svTestData"))
}

stats <- function (object, ...)
	UseMethod("stats")

stats.svTestData <- function (object, ...)
{
    if (!is.svTestData(object))
        stop("'object' must inherit from 'svTestData'")
    Stats <- attr(object, "stats")
    Table <- table(object$kind)
    ## Update the table with the total number of test
    Kinds <- c(Stats["tests"] - sum(Table[2:4], na.rm = TRUE), Table[2:4])
    names(Kinds) <- names(Table)
    ## Return a list with the table of kinds and the total timing
    return(list(kind = Kinds, timing = Stats["timing"]))
}

print.svTestData <- function (x, all = FALSE, header = TRUE, file = "",
append = FALSE, ...)
{
    ## If there is a context attribute, print info about the tests
    cat("", file = file, append = append)
    Context <- attr(x, "context")
    if (!is.null(Context)) {
        unitStr <- if (Context["unit"] == "") "" else
            paste(" (in ", basename(Context["unit"]), ")", sep = "")
        Stats <- stats(x)
        if (isTRUE(header))
			cat("\n== ", Context["test"], unitStr, .formatTime(Stats$timing,
				secDigits = 1), ": ", as.character(.kindMax(x$kind)), "\n",
				Context["msg"], "\n", sep = "", file = file, append = TRUE)
        cat(paste(c("//Pass:", "Fail:", "Errors:"), Stats$kind[1:3],
            collapse = " "), "//\n\n", sep = "", file = file, append = TRUE)
       ## Don't print tests that succeed if !all
        if (!isTRUE(all)) X <- x[x$kind != "OK", ] else X <- x
    } else X <- x
    ## Print info about each individual filtered test
    if (nrow(X) > 0) {
        Res <- ifelse(X$res == "", "", paste("\n", X$res, sep = ""))
        cat(paste("* ", X$msg, ": ", X$call, .formatTime(X$timing,
            secDigits = 3), " ... ", as.character(X$kind), Res, sep = "",
            collapse = "\n"), file = file, append = TRUE)
    }
    return(invisible(x))
}

summary.svTestData <- function (object, header = TRUE, file = "",
append = FALSE, ...)
{
    ## If there is a context attribute, print info about the tests
    cat("", file = file, append = append)
    Context <- attr(object, "context")
    if (!is.null(Context)) {
        unitStr <- if (Context["unit"] == "") "" else
            paste(" (in ", basename(Context["unit"]), ")", sep = "")
        Stats <- stats(object)
        if (isTRUE(header))
			cat("\n== ", Context["test"], unitStr, .formatTime(Stats$timing,
			    secDigits = 1), ": ", as.character(.kindMax(object$kind)), "\n",
			    Context["msg"], "\n", sep = "", file = file, append = TRUE)
        cat(paste(c("//Pass:", "Fail:", "Errors:"), Stats$kind[1:3],
            collapse = " "), "//\n\n", sep = "", file = file, append = TRUE)
    }
    ## List tests that failed
    Items <- rownames(object)
    Fail <- object$kind == "**FAILS**"
    if (any(Fail)) {
        cat("=== Failures\n", file = file, append = TRUE)
        cat(paste("[", Items[Fail], "] ", object$msg[Fail], ": ",
            object$call[Fail], collapse = "\n", sep = ""), "\n\n",
            sep = "", file = file, append = TRUE)
    }
    ## List tests that produce errors
    Err <- object$kind == "**ERROR**"
    if (any(Err)) {
        cat("=== Errors\n", file = file, append = TRUE)
        cat(paste("[", Items[Err], "] ", object$msg[Err], ": ",
            object$call[Err], collapse = "\n", sep = ""), "\n\n",
            sep = "", file = file, append = TRUE)
    }
}


protocol_junit.svTestData <- function (object, ...)
{
	if (!require(XML, quietly = TRUE))
		return(invisible(FALSE))

	toValidXmlString <- function (s) {
		s <- gsub("&", "&amp;", s)
		s <- gsub("<", "&lt;", s)
		s <- gsub(">", "&gt;", s)
		s <- gsub('"', "&quot;", s)
		s <- gsub("'", "&apos;", s)
		return(s)
	}

	basename <- function (s) sub(".*/", "", s)

	Context <- attr(object, "context")
	Stats <- attr(object, "stats")
	result <- xmlNode('testcase', attrs = c(
		'classname' = basename(Context[['unit']]),
        'name' = toValidXmlString(Context[['test']]),
        'time' = object$timing))
	kind <- as.numeric(.kindMax(object$kind))  # TODO: use accessor
	elementName <- c(NA, 'failure', 'error', NA)[kind]
	if (!is.na(elementName)) {
		failureNode <- xmlNode(elementName, attrs = c(
            'type' = elementName,
            'message' = toValidXmlString(object$res)))  # TODO: use accessor
		result <- addChildren(result, kids = list(failureNode))
	}
	if (kind == 4)
		result <- addChildren(result, kids = list(xmlNode('skipped')))

	return(result)
}
