`$.gui` <- function (gui, x, args)
{
	if (!exists(x, envir = gui)) return(NULL)
	res <- get(x, envir = gui)
	if (is.function(res)) {
		res <- function (...) get(x, envir = gui)(..., gui = gui)
		if (!missing(args))
			res <- do.call(res, args, envir = parent.frame())
	}
	return(res)
}

print.gui <- function (x, ...)
{
	## Pretty printing of 'gui' objects
	if (x$name == ".GUI") cat("The default SciViews GUI (.GUI)\n") else
		cat("A SciViews user interface named ", x$name, "\n", sep = "")
	cat("using widgets from: ", paste(guiWidgets(x), collapse = ", "),
		"\n", sep = "")
	if (!guiAsk(x)) cat("(it is currently inactivated - ask == FALSE)\n")
	if (!is.null(x[["call"]])) {
		cat("* Last call: ", deparse(x[["call"]]), "\n" , sep = "")
		cat("* Last widgets used: ", x$widgets, "\n", sep = "")
		cmt <- comment(x$status)
		if (is.null(cmt)) {
			if (is.null(x$status)) {
				cat("* Last status: ok\n")
			} else {
				cat("* Last status: ", x$status, "\n", sep = "")
			}
		} else cat("* Last status: ", x$status, " (", cmt, ")\n", sep = "")
		cat("* Last result:\n")
		print(x$res)
	}
	return(invisible(x))
}

is.gui <- function (x)
{
	## Does this object inherits from 'gui'?
	return(inherits(x, "gui"))
}
