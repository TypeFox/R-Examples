fileEdit <- function (..., title = files, editor = getOption("fileEditor"),
file.encoding = "", template = NULL, replace = FALSE, wait = FALSE)
{
    ## Rework files, title and template
	files <- c(...)
	lf <- length(files)
	if (!length(files)) {
		warning("You must provide at least one file path")
		return(invisible(FALSE))
	}	
	title <- rep(as.character(title), len = lf)
	if (length(template))
		template <- rep(as.character(template), len = lf)
	
	## If the file(s) do not exist or must be replaced,
	## create them (possibly from template)
	toReplace <- (isTRUE(as.logical(replace)) | !file.exists(files)) 
	if (length(toReplace)) {
		newFiles <- files[toReplace]
		if (length(template)) template <- template[toReplace]
		for (i in 1:length(newFiles)) {
			if (!length(template) || !nzchar(template[i])) {
				file.create(newFiles[i])
			} else if (file.exists(template[i])) {
				file.copy(template[i], newFiles[i], overwrite = TRUE,
					copy.mode = FALSE)
			} else { # Template file not found!
				warning("Template file '", template[i],
					'" not found, starting from an empty file')
				file.create(newFiles[i])
			}
		}    		
	}
	files <- normalizePath(files)
	
	## Manage file encoding
	if (nzchar(file.encoding) && file.encoding != "native.enc") {
        tfile <- file
        for (i in seq_along(file)) {
            tfile <- tempfile()
            con <- file(file[i], encoding = file.encoding)
            writeLines(readLines(con), tfile)
            close(con)
            file[i] <- tfile
        }
    }
	
	## There are a few shortcuts for editors that need to be expanded
	if (length(editor) && is.character(editor))
		editor <- switch(tolower(editor),
			textedit = "open -e -n -W \"%s\"",
			textwrangler = "edit --wait --resume \"%s\"",
			bbedit = "edit --wait --resume \"%s\"",
			editor)
	
	## Fallback to "editor", in case no fileEditor is provided
	if (!length(editor)) {
		editor <- getOption("editor")
	} else if (!isWin() && is.character(editor) && !grepl("%s", editor)) {
		cmd <- paste('which ', '"', editor, '"', sep = "")
		if (!length(system(cmd, intern = TRUE))) {
			## Fall back to the default editor (if any)
			editor <- getOption("editor")
		}	
	}
	
	## If not in interactive mode, or expressly no editor provided
	## We don't edit!
	if (!interactive() || !length(editor) ||
		(!is.function(editor) && !nzchar(editor))) {
		## Do nothing, issue, a warning!
		warning("Cannot edit files: no editor or not in interactive mode")
		return(invisible(FALSE))
	}
	
	## Special cases... where we prefer the internal editor
	## Note: just change editor a little bit to make sure to avoid internal!
	wait <- isTRUE(as.logical(wait))
	if (is.character(editor) &&
		editor %in% c("notepad", "internal", "vi", "open -e -n -W \"%s\"")) {
		done <- FALSE
		## 1) JGR
		if (isJGR()) {
			for (i in 1:lf)
				.fileEditJGR(files[i], title = title[i], wait = wait) 
			done <- TRUE
		## 2) Windows Rgui
		} else if (isRgui()) {
			for (i in 1:lf)
				.fileEditRgui(files[i], title = title[i], wait = wait)	
			done <- TRUE
		## 3) R.app and wait == FALSE (we cannot wait the end of edition using
		##    the internal R.app editor!) 
		} else if (isAqua() && !wait) {
			## Note that, here, the editor in use is the one defined in the
			## R.app preference dialog box!
			for (i in 1:lf)
				file.edit(files[i], title = title[i], fileEncoding = "")
			done <- TRUE
		}
		if (done) return(invisible(TRUE))
	}
	
	## In any other case, we use the defined editor
	if (is.function(editor)) {
		## Here, we need a special editor function that is able to wait!
		res <- try(editor(file = file, title = title, wait = wait),
			silent = TRUE)
	} else {
		## Construct the command...
		if (grepl("%s", editor)) {
			cmds <- sprintf(editor, files)
		} else {
			cmds <- paste('"', editor, '" "', files, '"', sep = "")
		}
		if (isMac()) msg <- "'... Close the editor (Cmd-Q) to continue!" else
			msg <- "'... Close the editor to continue!"
		for (i in 1:length(cmds)) {
			if (wait) message("Editing the file '", basename(files[i]), msg)
			flush.console()
			if (isWin()) {
				res <- try(system(cmds[i], ignore.stdout = TRUE,
					ignore.stderr = TRUE, wait = wait, minimized = FALSE,
					invisible = FALSE, show.output.on.console = FALSE),
					silent = TRUE)
			} else {
				res <- try(system(cmds[i], ignore.stdout = TRUE,
					ignore.stderr = TRUE, wait = wait), silent = TRUE)
			}
			if (inherits(res, "try-error")) break
		}
	}
	if (inherits(res, "try-error")) {
		warning(as.character(res))  # Transform the error into a warning
		return(invisible(FALSE))
	} else return(invisible(TRUE))
}

.fileEditJGR <- function (file, title = file, wait = FALSE) 
{
    ## Check arguments
	file <- as.character(file)
	if (length(title) != 1) stop("Only one item for 'file' is accepted")
	title <- as.character(title)
	if (length(title) != 1) stop("Only one item for 'title' is accepted")
	
	## Check if JGR is operating
	if (!isJGR()) {
        message(".fileEditJGR() cannot be used outside JGR.\n")
        return(invisible(NULL))
    }
	
	## Create a new editor window and open the file in it
	## Note that, if JGR is loaded, rJava is there too. So .jnew is available!
	editor <- rJava::.jnew('org/rosuda/JGR/editor/Editor', as.character(file)[1])
	## Set the title
	if (title != file) editor$setTitle(title)
	
	## Do we wait that the file is edited?
	if (isTRUE(as.logical(wait))) {
		message("Editing file '", basename(file),
			"'... Close the editor to continue!")
		while(editor$isVisible()) {
			editor$setState(0L)  # Make sure it is not iconized
			editor$toFront()     # Make the editor the frontmost window
			Sys.sleep(0.3)         # Wait for 0.3 sec
		}
	}
	invisible(editor)
} 

.fileEditRgui <- function (file, title = file, wait = FALSE) 
{
    ## Avoid errors in R CMD check about missing getWindowsHandles() function
	if (!isWin()) getWindowsHandles <- function(...) return()
	
	## Check arguments
	file <- as.character(file)
	if (length(title) != 1) stop("Only one item for 'file' is accepted")
	title <- as.character(title)
	if (length(title) != 1) stop("Only one item for 'title' is accepted")
	
	## Check if we are in RGui
	if (!isRgui()) {
        message(".fileEditRgui() cannot be used outside Rgui.\n")
        return(invisible(NULL))
    }
	
	## Edit file in an Rgui internal editor and track its existence
	## if wait == TRUE
	hdl <- getWindowsHandles()
	file.edit(file, title = title, editor = "internal", fileEncoding = "")
	hdl2 <- getWindowsHandles()
	editor <- hdl2[!hdl2 %in% hdl]
	
	## Do we wait that the file is edited?
	if (isTRUE(as.logical(wait)) && length(editor) == 1) {
		message("Editing file '", basename(file),
			"'... Close the editor to continue!")
		flush.console()
		while(editor %in% getWindowsHandles(minimized = TRUE))
			Sys.sleep(0.3)         # Wait for 0.3 sec
	}

	invisible(editor)
} 
