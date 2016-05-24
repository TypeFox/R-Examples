objInfo <- function (id = "default", envir = .GlobalEnv, object = "",
path = NULL)
{
	## Get a tooltip information for an object (for mouseover method)

	## Format envir as character (use only first item provided!)
	if (is.environment(envir)) envir <- deparse(substitute(envir))
	if (is.numeric(envir)) envir <- search()[envir[1]]
	envir <- as.character(envir)[1]

	## Possibly call a custom function .objInfo() in SciViews:TempEnv
	CmdFun <- getTemp(".objInfo", mode = "function")
    if (!is.null(CmdFun)) {  # We call a custom function
		Info <- CmdFun(id = id, envir = envir, object = object)
	} else if (object == "") {  # An environment...
		Info <- switch(envir,
			.GlobalEnv = paste(c("Global environment\n", capture.output(gc())),
				collapse = "\n"),
			`SciViews:TempEnv` = "SciViews temporary variables environment",
			RcmdrEnv = "R Commander temporary variables environment",
			`tools:RGUI` = "R.app tools environment",
			Autoloads = "R autoloading objects environment",
			if (regexpr("^package:", envir) > -1) {
				pkg <- sub("^package:", "", envir)
				paste(library(help = pkg, character.only = TRUE)$info[[1]],
					collapse = "\n")
			} else if (envir %in% search()) {
				paste("'", envir, "' environment", sep = "")
			} else ""
		)
	} else {  # An object...
		if (!exists(object, where = envir)) return(invisible(""))
		obj <- get(object, pos = envir)
		## The info is simply a str() representation of the object
		## We need to capture output
		Info <- capture.output(str(obj))
		## Add estimation of size for this object, if it is not a function
		if (!inherits(obj, "function")) {
			size <- object.size(obj)
			if (size > 1024*1024) {
				size <- paste("Estimated size:",
					format(size/1024/1024, digits = 3), "Mb")
			} else if (size > 1024) {
				size <- paste("Estimated size:",
					format(size/1024, digits = 3), "kb")
			} else size <- paste("Estimated size:",
					format(size, digits = 3), "bytes")
			Info[length(Info) + 1] <- size
		}
	}

	if (!is.null(path)) {
		## Save the data in a file
		if (path == "") path <- objDir()
		InfoFile <- file.path(path, paste("Info_", id, ".txt", sep = ""))
		cat(Info, collapse = "\n", file = InfoFile)
	}
	
	## Possibly call a .guiObjInfo function to pass the data to the GUI client
	CmdFun <- getTemp(".guiObjInfo", mode = "function")
    if (!is.null(CmdFun)) CmdFun(id = id, data = Info)

	## Return the info tooltip invisibly
	invisible(Info)
}
