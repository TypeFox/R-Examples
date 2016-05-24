addItems <- function (x, y, use.names = TRUE, replace = TRUE)
{
	if (replace) res <- c(y, x) else res <- c(x, y)
	if (use.names) {
		res <- res[!duplicated(names(res))]
	} else {
		res <- sort(unique(res))
	}
	res
}

addActions <- function (obj = ".svActions", text = NULL, code = NULL,
state = NULL, options = NULL, replace = TRUE)
{
	dat <- getTemp(obj, default = list())
	if (!inherits(dat, "list"))
		stop("'obj' should inherit from 'list'")

	## Make sure we return an svActions object
	class(dat) <- unique(c("svActions", class(dat)))

	## Add new actions characteristics to dat; make sure newdata are correct
	addData <- function(x, newdata, replace) {
		newnames <- names(newdata)
		if (is.null(newnames))
			stop("Data you add in actions must be a named character vector")
		newdata <- as.character(newdata)
		names(newdata) <- newnames
		x <- addItems(x, newdata, replace = replace)
		return(x)
	}
	if (!is.null(text)) dat$text <- addData(dat$text, text, replace)
	if (!is.null(code)) dat$code <- addData(dat$code, code, replace)
	if (!is.null(state)) dat$state <- addData(dat$state, state, replace)
	if (!is.null(options)) dat$options <- addData(dat$options, options, replace)

	## Reassign the modified values
	assignTemp(obj, dat)
	invisible(dat)
}

addIcons <- function (obj = ".svIcons", icons, replace = TRUE)
{
	## Get the list of icons
	icn <- getTemp(obj, default = character())
	if (!inherits(icn, "character"))
		stop("'obj' should inherit from 'character'")

	## Check that new icons are correctly formatted
	nicons <- names(icons)
	if (is.null(nicons))
		stop("Icons map you add must be a named character vector")
	icons <- as.character(icons)
	names(icons) <- nicons

	## Add new icons to it
	icn <- addItems(icn, icons, replace = replace)

	## Make sure we return an svIcons object
	class(icn) <- unique(c("svIcons", class(icn)))

	## Reassign the modified values
	assignTemp(obj, icn)
	invisible(icn)
}

addMethods <- function (methods)
{
	## Get the list of methods
	met <- getOption("svGUI.methods")
	if (!is.null(met)) methods <- addItems(met, methods, use.names = FALSE)
	options(svGUI.methods = sort(methods))
	invisible(methods)
}
