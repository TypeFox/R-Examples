batch <- function (items, fun, ...,	show.progress = !isAqua() && !isJGR(),
suppress.messages = show.progress, verbose = TRUE)
{
	if (!is.function(fun)) stop("'fun' must be a function")
	
	## Preparation of the batch process...
	owarn <- options(warn = 1) # Issue warnings immediatelly!
	on.exit(options(owarn))
	verbose <- isTRUE(as.logical(verbose))
	if (verbose) message("Running the batch process with ",
		deparse(substitute(fun)), "...")
	n <- length(items)
	if (n < 1) {
		warning("No items to process!")
		return(invisible(structure(FALSE, items = items, ok = logical(0))))
	}
	ok <- rep(NA, n) # A vector with results
	
	## Do we show progression?
	if (!isTRUE(as.logical(show.progress)))
		progress <- function (...) return() # Fake progress() function
	if (!isTRUE(as.logical(suppress.messages)))
		suppressMessages <- function (x) return(x) # Fake suppressMessages() fun
	
	## Run fun() for each item
	flush.console()
	for (i in 1:n) {
		progress(i, n)
		item <- items[i]
		ok[i] <- as.logical(suppressMessages(fun(item, ...)))[1]
		flush.console()
	}
	progress(n + 1, n) # Cancel progression message
	if (verbose) message("Processed successfully ", sum(ok, na.rm = TRUE),
		" items on ", n, " (see .last.batch)")
	## Record .last.batch variable in SciViews:TempEnv
	lastBatch <- structure(sum(ok, na.rm = TRUE) == n, items = items, ok = ok)
	assignTemp(".last.batch", lastBatch)
	invisible(lastBatch)
}
