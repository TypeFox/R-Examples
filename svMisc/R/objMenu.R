objMenu <- function(id = "default", envir = .GlobalEnv, objects = "",
sep = "\t", path = NULL)
{
	## TODO: look also in .required (now .Depends) in .GlobalEnv to determine
	## if one can detach a package
	## TODO: copy name to clipboard, send name to editor in menu

	## Get a context menu  for given object(s) or environment
	## It returns a matrix with the following columns:
	## - widget: "menu", "item", "sep" or "space"
	## - value: label of the menu or the item, or "-" for a separator
	## - tip: a short description of this menu entry
	## - code: the command to issue. Precede it with <<<h>>> to run the command
	##   silently and by <<<H>>> to run it silently and disconnect immediately
	##   Use a substitution placeholder for replacement in value, tip and code.
	##   For instance, <<<obj>>>, <<<objname>>>, <<<pos>>>, <<<envir>>> and
	##   add an argument returning the string to place there in the call to
	##   .addStripbar()
	## - icon: the URL of the icon to use in the menu
	## - checked: if the menu entry should present a check mark
	## - disabled: if the menu entry is disabled
	## - hidden: if the menu entry is hidden
	## - options: further options that the GUI client can interpret
	## Notes:
	## - Selecting environments and objects inside environments is not allowed
	## - Selecting multiples environments is not supported
	## - Selectiong items in different environments is not allowed

	popup <- .createStripbar("popupbar")

	## Format envir as character (use only first item provided!)
	if (is.environment(envir)) envir = deparse(substitute(envir))
	if (is.numeric(envir)) envir <- search()[envir[1]]
	envir <- as.character(envir)[1]
	## Get the current position in the search path for envir
	pos <- match(envir, search(), nomatch = -1)
	if (pos < 1) return(invisible(popup))  # Not found!

	## Do we replace envir or not?
	if (envir == ".GlobalEnv") Envir <- "<<<envir>>>" else Envir <- envir

	## Look at 'objects'
	if (length(objects) > 1) {  # We have a multiple selection
		objType <- "multiple"
		## Add 'save' and 'remove' menu entries
		objsPar <- paste(objects, collapse = ", ")
		popup <- .addStripbar(popup, c("save", "remove"), obj = objsPar,
			envir = Envir)

	} else {  # Only one item is selected
		if (objects == "") {  # This is a menu for the environment
			## The menu is different depending if we are in .GlobalEnv or not
			if (envir == ".GlobalEnv") {
				objType <- ".GlobalEnv"
				objPar <- NULL  # Nothing special to propose
				popup <- .addStripbar(popup, c(
						"load",
						"source",
						"import",
						"sep_",
						"detach (d)"
					), envir = envir)

			} else if (regexpr("^package:", envir) > -1) {  # A package envir
				objType <- "package"
				objPar <- sub("^package:", "", envir)
				req <- getOption("required.packages")
				if (is.null(req))  # Use list of known sensible environments
					req <- c("package:base", "package:methods", "package:datasets",
						"package:utils", "package:grDevices", "package:graphics",
						"package:stats", "package:tcltk", "package:svMisc",
						"package:svGUI", "package:svSocket", "package:svViews",
						"package:svIO", "package:svIDE", "package:svDialogs",
						"package:svWidgets", "package:tcltk2")
				## Make sure that "package:base" is in the list
				req <- c("package:base", req)
				detachable <- !(envir %in% req)
				if (detachable) state <- "" else state <- " (d)"
				popup <- .addStripbar(popup, c(
						"pkgInfo",
						"sep_",
						paste("detach", state, sep = ""),
						paste("detachUnload", state, sep = "")
					), package = objPar, envir = Envir)

			} else {  # Another environment than .GlobalEnv, but not a package
				objType <- "environment"
				req <- getOption("required.environments")
				if (is.null(req))  # Use list of known sensible environments
					req <- c("SciViews:TempEnv", "RcmdrEnv", "Autoloads", "tools:RGUI")
				detachable <- !(envir %in% req)
				if (detachable) state <- "" else state <- " (d)"
				objPar <- detachable
				popup <- .addStripbar(popup, c(paste("detach", state, sep = "")),
					envir = Envir)
			}

		} else {  # This is a menu for an object
			if (!exists(objects, where = envir)) return(invisible(popup))
			objType <- "object"
			obj <- get(objects, pos = envir)
			objPar <- class(obj)

			## Look if there is an help file and examples for this object
			hlp <- isHelp(objects)
			names(hlp) <- NULL
			## Then we add 'help' and 'example' at the top of the menu
			popup <- .addStripbar(popup, c(
				ifelse(hlp, c("help", "example"),
				c("help (d)", "example (d)")), "sep_"), obj = objects)

			## Add a series of standard functions/methods for this object
			mets <- listMethods(class = objPar)
			popup <- .addStripbar(popup, c(
					"Functions_.",
					"._print",
					ifelse("summary" %in% mets, "._generic", "._generic (d)")
				), obj = objects, fun = "summary")
			popup <- .addStripbar(popup, c(
					ifelse("plot" %in% mets, "._generic", "._generic (d)"),
					ifelse(!is.null(names(obj)), "._names", "._names (d)"),
					"._str"
				), obj = objects, fun = "plot")
			## Eliminate 'print', 'show', 'summary' and 'plot' from mets
			mets <- mets[!(mets %in% c("print", "show", "summary", "plot"))]
			## Add the other methods from mets in the menu
			if (length(mets) > 0) {
				popup <- .addStripbar(popup, "._sep_")
				for (met in mets)
					popup <- .addStripbar(popup, "._generic", obj = objects,
						fun = met)
			}

			## Add a &View menu and at least one entry: &Default view
			enabledViews <- "svViews" %in% search()	# We need this package
			if (enabledViews) {
				popup <- .addStripbar(popup, c(
						"View_.",
						"._viewDef"
					), obj = objects)
				## Complete the '&View' submenu with more entries
				viewTypes <- listTypes("view", objPar)
				for (viewType in viewTypes)
					popup <- .addStripbar(popup, "._view", obj = objects,
						type = viewType)
				popup <- .addStripbar(popup, "report", obj = objects)
			} else {  # Give the opportunity to load the svViews package
				popup <- .addStripbar(popup, c(
						"View_.",
						"._viewDef (d)",
						"._require",
						"report (d)"
					), obj = objects, pkg = "svViews")
			}

			## Add &edit/fix and &save
			popup <- .addStripbar(popup, c(
					ifelse(envir == ".GlobalEnv", "edit", "fix"),
					"save"
				), obj = objects, envir = Envir)

			## If we are in .GlobalEnv and the object is a data.frame or a list,
			## then I can attach/detach it
			if (envir == ".GlobalEnv" && inherits(obj, "list")) {
				## Is this object already attached (present in the search path)
				if (objects %in% search()) {
					popup <- .addStripbar(popup, c("detach", "reattach", "sep_"),
						obj = objects)
				} else {  # The object is not attached yet...
					popup <- .addStripbar(popup, c("attach", "sep_"),
						obj = objects)
				}
			}

			## Create the '&copy' submenu with all possible entries
			## plus &Export and &Remove
			enabledCopy <- "svIO" %in% search()	# We need this package
			if (enabledCopy) {
				popup <- .addStripbar(popup, c(
					"Copy_.",
					"._copyDef"
					), obj = objects)
				## Complete the '&Copy' submenu with more entries
				copyTypes <- listTypes("copy", objPar)
				for (copyType in copyTypes)
					popup <- .addStripbar(popup, "._copy", obj = objects,
						type = copyType)
					popup <- .addStripbar(popup, c(
						"export",
						"sep_",
						"remove"
					), obj = objects, envir = Envir)
			} else {  # Give the opportunity to load the svIO package
					popup <- .addStripbar(popup, c(
						"Copy_.",
						"._copyDef (d)",
						"._require",
						"export (d)",
						"sep_",
						"remove"
					), obj = objects, envir = Envir, pkg = "svIO")
			}
		}
	}  # Done (menu construction)

	## Possibly call a custom function .objMenu() in SciViews:TempEnv
	## to finalize the menu
	## Depending on objType, we send something we already calculated to avoid
	## doing it twice:
	## objType   -> objPar
	## .GlobalEnv   NULL
	## package      Package name
	## environment  is it detachable or not?
	## object       class of the object
	## multiple     string with the multiple objects separated by ', '
	CmdFun <- getTemp(".objMenu", mode = "function")
    if (!is.null(CmdFun))
		popup <- CmdFun(popup, id = id, envir = envir, pos = pos,
			objects = objects, path = path, objType = objType, objPar = objPar)

	if (!is.null(path)) {
		## Save the data in a file
		if (path == "") path <- objDir()
		MenuFile <- file.path(path, paste("Menu_", id, ".txt", sep = ""))
		write.table(popup, file = MenuFile, sep = sep, row.names = FALSE,
			col.names = FALSE)
	}

	## TODO: make this more flexible
	## Possibly call a .guiObjMenu function to pass the data to the GUI client
	CmdFun <- getTemp(".guiObjMenu", mode = "function")
    if (!is.null(CmdFun)) CmdFun(id = id, data = popup)

	## Return the menu specification invisibly
	invisible(popup)
}
