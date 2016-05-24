objBrowse <- function (id = "default", envir = .GlobalEnv, all.names = NULL,
pattern = NULL, group = NULL, sep = "\t", path = NULL, regenerate = FALSE)
{
    ## Maintain files for remote Object Browser
    ## If four first parameters are NULL, use cached version of these parameters,
	## or default values

	## Format envir as character (use only first item provided!)
	if (is.environment(envir)) envir <- deparse(substitute(envir))
	if (is.numeric(envir)) envir <- search()[envir[1]]
	envir <- as.character(envir)[1]
	## Get the current position in the search path for envir
	pos <- match(envir, search(), nomatch = -1)
	if (pos < 1) {
		pos <- 1  # NOT FOUND, use .GlobalEnv
		envir = ".GlobalEnv"
	}

    if (!is.null(path)) {
		## Does the directory exists?
		if (path == "") path <- objDir()
		if (!file.exists(path) || !file.info(path)$isdir) {
			#unlink(path)
			if (!dir.create(path))
				stop("Impossible to create the Object Browser 'path' directory!")
		}
	}

	## Control that 'Search.txt' is up-to-date
	ChangedSearch <- objSearch(path = path, compare = !regenerate)

	## Make sure id is character
	id <- as.character(id)[1]
	if (id == "") id <- "default"

	## Get the five parameters pos, envir, all.names, pattern & group
    allPars <- getTemp(".guiObjParsCache", default = NULL)
    if (!is.null(allPars)) {
    	Pars <- allPars[[id]]
    } else {
    	Pars <- list(pos = 1, envir = ".GlobalEnv", all.names = FALSE,
			pattern = "", group = "")
    	assignTemp(".guiObjParsCache", list())	# Create the list
    }
    if (is.null(Pars))
		Pars <- list(pos = 1, envir = ".GlobalEnv", all.names = FALSE,
			pattern = "", group = "")
    ## Possibly change some parameters (and make sure they are valid!)
	ParsChanged <- FALSE
	if (!is.null(pos)) {
        ParsChanged <- TRUE
		Pars$pos <- as.integer(pos[1])
		Pars$envir <- envir
		if (Pars$pos < 1) {
			ParsChanged <- TRUE
			Pars$pos <- 1
			Pars$envir <- ".GlobalEnv"
		}
	}
	if (Pars$pos > length(search())) {
		ParsChanged <- TRUE
		Pars$pos <- 1
		Pars$envir <- ".GlobalEnv"
	}
	## Track possible changes in the search path
	if (is.na(match(Pars$envir, search()))) {
		ParsChanged <- TRUE
		Pars$pos <- 1
		Pars$envir <- ".GlobalEnv"
	}
	if (match(Pars$envir, search()) != Pars$pos) {
		ParsChanged <- TRUE
		Pars$pos <- match(Pars$envir, search())
	}
	## Track changes in the options
	if (!is.null(all.names)) {
		ParsChanged <- TRUE
		Pars$all.names <- as.logical(all.names[1])
	}
	if (!is.null(pattern)) {
		ParsChanged <- TRUE
		Pars$pattern <- as.character(pattern[1])
	}
	if (!is.null(group)) {
		ParsChanged <- TRUE
		Pars$group <- as.character(group[1])
	}
	## Write a cached version of these parameters in SciViews:TempEnv
	allPars <- getTemp(".guiObjParsCache", default = list())
	allPars[[id]] <- Pars
	assignTemp(".guiObjParsCache", allPars)

    ## Control that 'List_<id>.txt' is up-to-date, but only if pos == 1 or
	## envir is not a package or regenerate or Pars or Search have changed
	## to limit the work done on workspaces that are unlikely to have change
	if (ParsChanged || regenerate || (ChangedSearch != "") || (Pars$pos == 1) ||
		(regexpr("^package:", envir) == -1)) {
		ChangedList <- objList(id = id, envir = Pars$pos,
			all.names = Pars$all.names, pattern = Pars$pattern,
			group = Pars$group, path = path, compare = !regenerate,
			sep = sep)

		ChangedList <- if(!is.null(nrow(ChangedList)) && nrow(ChangedList) > 0) {
			apply(ChangedList, 1, paste, collapse = sep)
		} else ""

	} else ChangedList <- ""

	## We return the data, or indication that the data have changed to the client
	Data <- ""
	if (length(ChangedSearch) > 1 || ChangedSearch != "") {
		Data <- "<<<search>>>\n"
		if (is.null(path)) Data <- paste(Data,
			paste(ChangedSearch, collapse = sep), "\n", sep = "")
	}
	if (length(ChangedList) > 1 || ChangedList != "") {
		Data <- paste(Data, "<<<list>>>", sep = "")
		if (is.null(path)) Data <- paste(c(Data, ChangedList),
			collapse = "\n")
	}

	## Possibly call a .guiObjBrowse function to pass the data to the GUI client
	CmdFun <- getTemp(".guiObjBrowse", mode = "function")
    if (!is.null(CmdFun)) CmdFun(id = id, data = Data)
	## Return the data invisibly
	invisible(Data)
}
