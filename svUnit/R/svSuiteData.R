is.svSuiteData <- function (x)
{
	## It this a svSuiteData object
	return(inherits(x, "svSuiteData"))
}

stats.svSuiteData <- function (object, ...)
{
    if (!is.svSuiteData(object))
        stop("'object' must inherit from 'svSuiteData'")
    ## Get the list of tests
    Tests <- ls(object)
    if (length(Tests) == 0) {
        ## The object is empty!
        Res <- data.frame(kind = .kind(logical()), timing = numeric(),
            time = numeric(), unit = character(), tag = character(),
            msg = character(), stringsAsFactors = FALSE)
    } else {
        ## Functions to get data for each test
        getKind <- function(x) .kindMax(x$kind)
        getTiming <- function(x) attr(x, "stats")["timing"]
        getTime <- function(x) attr(x, "time")
        getContext <- function(x, item) attr(x, "context")[[item]]
        Res <- data.frame(
            kind = rev(sapply(object, getKind)),
            timing = rev(sapply(object, getTiming)),
            time = structure(rev(sapply(object, getTime)),
                class = c("POSIXt", "POSIXct")),
            unit = rev(sapply(object, getContext, "unit")),
            msg = rev(sapply(object, getContext, "msg")),
            stringsAsFactors = FALSE)
    }
    return(Res)
}

metadata <- function (object, ...)
	UseMethod("metadata")

metadata.svSuiteData <- function (object,
fields = c("R.version", "sessionInfo", "time", "description"), ...)
{
    ## Extract metadata information from a 'svSuiteData' object
	if (!is.svSuiteData(object))
		stop("'object' must inherit from 'svSuiteData'")
	## Return a list with all metadata elements found
	fields <- paste(".", fields, sep = "")
    Res <- list()
	for (F in fields)
		Res[[F]] <- object[[F]]
	return(Res)
}

print.svSuiteData <- function (x, all = FALSE, file = "", append = FALSE, ...)
{
    if (!is.svSuiteData(x))
        stop("'x' must inherit from 'svSuiteData'")
    Tests <- ls(x)
    if (length(Tests) == 0) {
        cat("No test records!\n", file = file, append = append)
    } else {
        ## Print general information about the tests
        Stats <- stats(x)
		Tests <- rownames(Stats)	# To make sure we use the same!
        Timing <- .formatTime(sum(Stats$timing, na.rm = TRUE), secDigits = 1)
		cat("= A svUnit test suite", Timing, " with:\n\n", sep = "",
			file = file, append = append)
		cat(paste("* ", Tests, " ... ", as.character(Stats$kind), "",
			sep = "", collapse = "\n"),
            "\n\n", sep = "", file = file, append = TRUE)

        ## Print detailed information about each test
        for (Test in Tests)
            print(x[[Test]], all = all, file = file, append = TRUE, ...)
    }
    return(invisible(x))
}

summary.svSuiteData <- function (object, ...)
    protocol_text.svSuiteData(object, ...)

protocol <- function (object, type = "text", file = "", append = FALSE, ...)
	UseMethod("protocol")

protocol.default <- function (object, type = "text", file = "", append = FALSE, ...)
	get(paste("protocol", type[1], sep = "_"))(object, file = file, append = append, ...)

protocol.svSuiteData <- function (object, type = "text", file = "", append = FALSE, ...)
	get(paste("protocol", type[1], sep = "_"))(object, file = file, append = append, ...)

protocol_text <- function (object, file = "", append = FALSE, ...)
	UseMethod("protocol_text")

protocol_text.svSuiteData <- function (object, file = "", append = FALSE, ...)
{
    if (!is.svSuiteData(object))
        stop("'object' must inherit from 'svSuiteData'")
    Tests <- sort(ls(object))
    if (length(Tests) == 0) {
        cat("No test records!\n", file = file, append = append)
    } else {
        ## Print general information about the tests
        Stats <- stats(object)
		Tests <- rownames(Stats)	# To make sure we use the same!
		Timing <- .formatTime(sum(Stats$timing, na.rm = TRUE), secDigits = 1)
		cat("= A svUnit test suite", Timing, " with:\n\n", sep = "",
			file = file, append = append)
        cat(paste("* ", Tests, " ... ", as.character(Stats$kind), "",
			sep = "", collapse = "\n"),
            "\n\n", sep = "", file = file, append = TRUE)

        ## Summarize each test
        for (Test in Tests)
            summary(object[[Test]], file = file, append = TRUE)
    }
}

protocol_junit <- function (object, ...)
	UseMethod("protocol_junit")

protocol_junit.svSuiteData <- function (object, file = "", append = FALSE, ...)
{
	if (!is.svSuiteData(object))
		stop("'object' must inherit from 'svSuiteData'")
	if(!require(XML, quietly = TRUE))
		return(invisible(FALSE))

	Tests <- sort(ls(object))
	if (length(Tests) > 0 && inherits(object[[Tests[1]]], "svSuiteData")) {
		## This is a set of suites (containing svSuiteData)
		root <- xmlNode('testsuites')
	} else {
		## This is a single suite (containing svTestData)
		root <- xmlNode('testsuite')
	}

	with(stats(object), addAttributes(root, name = NULL, tests = length(Tests),
		errors = sum(kind == '**ERROR**'), failures = sum(kind == '**FAILS**'),
        skip = sum(kind == 'DEACTIVATED')))
	for (Test in Tests)
		root <- addChildren(root,
			kids = list(protocol_junit(object[[Test]], append = TRUE)))

	## Decide whether to return the xml node or write the xml file
	if (isTRUE(append)) {
		return(root)
	} else {
		saveXML(root, file)
		return(invisible(TRUE))
	}
}
