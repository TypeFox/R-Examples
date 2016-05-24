.onLoad <- function (lib, pkg) {
	.initialize()

	## Determine where to find the preferred file editor for fileEdit()
	if (is.null(getOption("fileEditor"))) {
		if (interactive()) {
			if (isWin()) {
				## First check if Notepad++ is installed in default location...
				pfdir <- Sys.getenv("ProgramFiles")
				if (pfdir == "") pfdir <- "c:\\program files"
				fEdit <- paste(pfdir, "\\Notepad++\\notepad++.exe", sep = "")
				## ... otherwise, fallback to notepad.exe
				if (!file.exists(fEdit)) fEdit <- "notepad"  # Default for Rterm
			} else if (isMac()) {
				## Note that, in R.app, one cannot edit and wait for it is done
				## So, I always define a different editor, but fall back to
				## internal R.app if internal.if.possible is TRUE
				## First check if 'edit' is there
				## (open files in TextWrangler or BBEdit)
				if (length(suppressWarnings(system("which edit",
					intern = TRUE)))) {
					fEdit <- "textwrangler"
				} else { # Fall back to the default text editor
					## Note: use "open -e -n -W \"%s\"" to force use of TextEdit
					fEdit <- "textedit"
				}			
			} else { # This is unix
				## First check if gedit or kate is there... This is different
				## from file.edit() and editor that looks directly to EDITOR!
				if (length(suppressWarnings(system("which gedit",
					intern = TRUE)))) {
					fEdit <- "gedit"
				} else if (length(suppressWarnings(system("which kate",
					intern = TRUE)))) {
					fEdit <- "kate"
				} else { # Fall back to something that is likely to be installed
					## Look id EDITOR or VISUAL environment variable is defined
					fEdit <- Sys.getenv("EDITOR")
					if (nzchar(fEdit)) fEdit <- Sys.getenv("VISUAL")
					if (nzchar(fEdit)) fEdit <- "vi"
				}
			}
			options(fileEditor = fEdit)
		} else options(fileEditor = "") # Inactivate it!
	}	
}

.initialize <- function (replace = TRUE)
{
	## Create .svActions if it does not exists yet
	.svActions <- list()
	class(.svActions) <- unique(c("svActions", class(.svActions)))
	assignTemp(".svActions", .svActions, replace.existing = FALSE)

	## Define actions we need for the object browser menus
	addTemp(".svActions", "text", c(
		load =      gettext("Load...\nLoad R objects"),
		source =    gettext("Source...\nSource R code"),
		save =      gettext("Save as...\nSave to a file"),
		import =    gettext("Import...\nImport data in R"),
		export =    gettext("Export...\nExport data to a file"),
		report =    gettext("Report...\nPrepare a report for this object"),
		setwd =     gettext("Set Working dir...\nChange current R working directory"),
		print =     gettext("Print or show\nPrint or show the content of the object"),
		generic =   gettext("<<<fun>>>()\nApply method <<<fun>>>() to the object"),
		names =     gettext("Names\nNames of variables contained in the object"),
		str =       gettext("Str\nCompact str() representation of an object"),
		help =      gettext("Help\nHelp on an object"),
		example =   gettext("Example\nRun examples for this object"),
		edit =      gettext("Edit\nEdit an object"),
		fix =       gettext("Fix\nFix an R object"),
		pkg =       gettext("Load package(s)\nLoad one or several R packages"),
		remove =    gettext("Remove\nRemove (permanently!) one or several objects from memory"),
		require =   gettext("Require <<<pkg>>>\nRequire the package <<<pkg>>>"),
		attach =    gettext("Attach\nAttach an object to the search path"),
		detach =    gettext("Detach\nDetach an object or package from the search path"),
		detachUnload = gettext("Detach and unload\nDetach a package from the search path and unload it"),
		reattach =  gettext("Reattach\nReattach an object to the search path"),
		pkgInfo =   gettext("Package info\nShow detailed information for this package"),
		viewDef =   gettext("View (default)\nDefault view for this object"),
		view =      gettext("View <<<type>>>\nDisplay a '<<<type>>>' view for this object"),
		copyDef =   gettext("Copy (default)\nCopy this object to the clipboard (default format)"),
		copy =      gettext("Copy <<<type>>>\nCopy this object to the clipboard in '<<<type>>>' format"),
		Functions = gettext("Functions\nGeneric functions and methods"),
		View =      gettext("View\nView the object"),
		Copy =      gettext("Copy\nCopy the object to the clipboard")
	), replace = replace)

	addTemp(".svActions", "code", c(
		load =     "guiLoad([[[pos = \"<<<envir>>>\"]]])",
		source =   "guiSource([[[pos = \"<<<envir>>>\"]]])",
		save =     "guiSave(<<<obj>>>[[[, pos = \"<<<envir>>>\"]]])",
		import =   "guiImport()",
		export =   "guiExport(<<<obj>>>)",
		report =   "guiReport(<<<obj>>>)",
		setwd =    "guiSetwd([[[<<<dir>>>]]])",
		print =    "<<<obj>>>",
		generic =  "[[[<<<var>>>> <- ]]]<<<fun>>>(<<<obj>>>)",
		names =    "names(<<<obj>>>)",
		str =      "str(<<<obj>>>)",
		help =     "help(<<<obj>>>)",
		example =  "example(<<<obj>>>)",
		edit =     "<<<obj>>> <- edit(<<<obj>>>)",
		fix =      "fix(<<<obj>>>)",  # There is no guarantee we fix the right one!
		pkg =      "[[[<<<res>>> <- ]]]pkg(\"<<<pkgs>>>\")",
		remove =   "rm(<<<obj>>>[[[, pos = \"<<<envir>>>\"]]])",
		require =  "[[[<<<res>>> <- ]]]require(<<<pkg>>>)",
		attach =   "attach(<<<obj>>>)",
		detach =   "detach(<<<envir>>>)",
		detachunload = "detach(<<<envir>>>, unload = TRUE)",
		reattach = "detach(<<<obj>>>); attach(<<<obj>>>)",
		pkgInfo =  "<<<H>>>library(help = <<<package>>>)",
		viewDef =  "view(<<<obj>>>)",
		view =	   "view(<<<obj>>>, type = \"<<<type>>>\")",
		copyDef =  "copy(<<<obj>>>)",
		copy =	   "copy(<<<obj>>>, type = \"<<<type>>>\")"
	), replace = replace)

	addTemp(".svActions", "state", c(
		viewDef = "d",
		copyDef = "d"
	), replace = replace)

	addTemp(".svActions", "options", c(
		generic = ""
	), replace = replace)

	## If the option svGUI.methods is not defined, give reasonable default values
	## Those are methods that can be applied to many objects without providing
	## additional argument and that will be added automatically to objects'
	## context menu in GUIs (print and show are not included, because we know
	## they must exist for all objects)
	## Rem: use addMethods() if you just want to add methods to this list
	if (is.null(getOption("svGUI.methods")))
		options(svGUI.methods = c("AIC", "anova", "confint", "BIC", "formula",
			"head", "hist", "logLik", "plot", "predict", "residuals", "summary",
			"tail", "vcov"))
}

.createStripbar <- function (type = c("menubar", "popupbar", "toolbar", "buttonbar", "statusbar"))
{
	type <- match.arg(type)
	strp <- data.frame(widget = character(), value = character(),
		tip = character(), code = character(), icon = character(),
		checked = logical(), disabled = logical(), hidden = logical(),
		options = character(), stringsAsFactors = FALSE)
	class(strp) <- unique(c("svStripbar", "svStrip", class(strp)))
	attr(strp, "type") <- type
	return(strp)
}

.addStripbar <- function (strip, widgets, gui = getOption("svGUI.name"), actions = NULL, ...,
	icons = NULL)
{
	if (!inherits(strip, "svStripbar"))
		stop("'strip' must be a 'svStripbar' object")
	Type <- attr(strip, "type")
	if (is.null(Type)) Type <- "menubar"	# Default value
	## Extract possible state information from widgets
	pos <- regexpr(" *[(][cCuUdDeEhHvV]+[)] *$", widgets)
	pos[pos == -1] <- 1000000
	state <- substr(widgets, pos, 1000000)
	## Clean up state to keep only digits
	state <- tolower(sub("^ *[(]([cCuUdDeEhHvV]+)[)] *$", "\\1", state))
	widgets <- substr(widgets, 1, pos - 1)

	## If widgets is a named character vector, just check it is
	## 'menu', 'item', 'sep' or 'space', otherwise, compile name = widget
	wnames <- names(widgets)
	if (is.null(wnames)) {
		wnames <- widgets
		## Determine widgets according to wnames
		widgets <- rep("item", length.out = length(wnames))
		widgets[regexpr("_[.]$", wnames) > -1] <- "menu"
		widgets[regexpr("_$", wnames) > -1] <- "sep"
		widgets[regexpr("__$", wnames) > -1] <- "space"
		names(widgets) <- wnames
	} else {
		## Name is provided => check content for 'menu', 'item', 'sep' or 'space'
		if (!all(widgets %in% c("menu", "item", "sep", "space")))
			stop("'widgets' must be 'menu', 'item', 'sep' or 'space'")
	}

	## Get tree hierarchy of the menus being the number of dots before '_'
	tree <- sub("^([.]+)_.*$", "\\1", wnames)
	tree[regexpr("^[.]+$", tree) == -1] <- ""
	tree <- gsub("[.]", "|", tree)

	## Clean up widget names
	wnames <- sub("^[.]+_", "", wnames)
	wnames <- sub("_[.|_]{0,1}$", "", wnames)

	## Guess some of the values from the names
	valuedef <- wnames
	valuedef[widgets == "sep"] <- "-"
	valuedef[widgets == "space"] <- "<->"

	## Collect together 'text', 'code', 'state' and 'options' from actions,
	## .svActions.[gui] and .svActions
	if (is.null(gui)) gui <- "___"	# Default value if no gui exists
	act.gui <- getTemp(paste(".svActions", gui, sep = "."),
		default = structure(list(), class = c("svActions", "list")))
	act <- getTemp(".svActions",
		default = structure(list(), class = c("svActions", "list")))
	## Collect together items
	deftext <- c(actions$text, act.gui$text, act$text)
	defcode <- c(actions$code, act.gui$code, act$code)
	defstate <- c(actions$state, act.gui$state, act$state)
	defoptions <- c(actions$options, act.gui$options, act$options)

	## Do the same for icons
	deficons <- c(icons,
		getTemp(paste(".svIcons", gui, sep = "."), default = character()),
		getTemp(".svIcons", default = character()))

	## The function used to replace placeholders in text and code
	replace <- function (x, ...) {
		## Do replacement for ... arguments
		args <- list(...)
		largs <- length(args)
		if (length(args) > 0) {
			nargs <- names(args)
			for (i in 1:largs)
				if (nargs[i] != "")
					x <- gsub(paste("<<<", nargs[i], ">>>", sep = ""),
						args[[i]], x)

		}
		## Eliminate optional parts where no replacement occured
		x <- gsub("\\[\\[\\[.*<<<.*>>>.*\\]\\]\\]", "", x)
		## Eliminate triple square brackets for optional parts we keep
		x <- gsub("\\[\\[\\[", "", x)
		x <- gsub("\\]\\]\\]", "", x)
		return(x)
	}

	## Compute the different elements we need
	## - text => value (first line) and tip (the rest)
	text <- replace(deftext[wnames], ...)
	text[is.na(text)] <- ""	# TODO: a better guess for menus, items and sep/space
	names(text) <- wnames
	pos <- regexpr("\n", text)
	pos[pos == -1] <- 1000000
	value <- substr(text, 1, pos - 1)
	value[value == ""] <- valuedef[value == ""]
	tip <- substr(text, pos + 1, 1000000)

	## Indicate menu hierarchy in value
	value <- paste(tree, value, sep = "")

	## - code
	code <- replace(defcode[wnames], ...)
	code[is.na(code)] <- ""
	names(code) <- wnames

	## - icon
	icon <- deficons[wnames]
	icon[is.na(icon)] <- ""
	names(icon) <- wnames

	## - options
	options <- defoptions[wnames]
	options[is.na(options)] <- ""
	names(options) <- wnames

	## Compute state => checked, disabled and hidden
	state2 <- defstate[wnames]
	state2[is.na(state2)] <- ""
	state <- paste(state, tolower(state2))
	checked <- ifelse(regexpr("c", state) > -1, TRUE, FALSE)
	disabled <- ifelse(regexpr("d", state) > -1, TRUE, FALSE)
	hidden <- ifelse(regexpr("h", state) > -1, TRUE, FALSE)

	## Create the data frame containing the data
	addstrip <- data.frame(widget = widgets, value = value, tip = tip,
		code = code, icon = icon, checked = checked, disabled = disabled,
		hidden = hidden, options = options, stringsAsFactors = FALSE)
	snames <- rownames(strip)
	## Add it to strip and return it
	strip <- rbind(strip, addstrip)
	rownames(strip) <- make.names(c(snames, wnames), unique = TRUE)
	## Make sure class and type are kept
	class(strip) <- unique(c("svStripbar", "svStrip", class(strip)))
	attr(strip, "type") <- Type
	return(strip)
}

#test <- c(
#	"load",
#	"sep_",
#	"File_.",
#	"._import (cdh)",
#	"._sep_",
#	"._export (h) ",
#	"._View_.",
#	".._viewDef (ch)",
#	".._view",
#	"._report",
#	"attach  (d)  "
#)
#pop <- .createStripbar("popupbar")
#.addStripbar(pop, test, obj = "testobj", type = "mytype")

## gettext() and hence gettextf() cannot retrieve messages ending with space
## in the "R" domain, because these functions stripe them out!
## This is a hack using ngettext() that uses unmodified version of the message
## Restriction: on the contrary to gettext(), .gettext() can translate only
## one message at a time, and default domain is changed to "R"
.gettext <- function (msg, domain = "R")
    ngettext(1, msg, "", domain = domain)

.gettextf <- function (fmt, ..., domain = "R")
    sprintf(ngettext(1, fmt, "", domain = domain), ...)
