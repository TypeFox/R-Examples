tkWinAdd <- function (name = "win1", parent = .TkRoot, title = NULL, pos = NULL,
bind.delete = TRUE, ...)
{
	## Verify the name is valid and not yet in use
	if (is.null(name) || !inherits(name, "character"))
		stop("'name' must be a character string!")
	name <- name[1]
	.guiWins <- getTemp(".guiWins")
	if (!is.null(.guiWins[[name]]))
		stop("name ", name, " is already used!")
	## It is fine. We can create the tktoplevel window
	win <- tktoplevel(parent = parent, ...)
	if (!inherits(win, "tkwin")) {
		tkdestroy(win)
		stop("A problem occured while creating the window!")
	} else class(win) <- c("tkguiWin", "guiWin", "gui", class(win))
	if (!is.null(title))
		tktitle(win) <- title[1]
	## Possibly position the window
	if (!is.null(pos))
		tkwm.geometry(win, pos)
    if (isTRUE(bind.delete)) {
		## Define action when clicking on 'X'
		##tkbind(win "<Destroy>", function() {tkWinDel(name); tkdestroy(win)})
		tkwm.protocol(win, "WM_DELETE_WINDOW", function() tkWinDel(name))
	}
	## Possibly change other characteristics of the window here...
	## ...

	## Record the window in the list
	.guiWins[[name]] <- win
	assignTemp(".guiWins", .guiWins)
	win
}

tkWinDel <- function (window)
{
	## Same action as tkdestroy(), but cares about deleting relating resources
	## (i.e., in .guiXXX in SciViews:TempEnv)
    win <- winGet(window)
	if (is.null(win))
		return(invisible(FALSE))	# window does not exist
    ## Delete it from the windows list
	.guiWins <- getTemp(".guiWins")
	.guiWins[[window]] <- NULL
	## Eliminate all related menus
	tkMenuDel(paste("$Tk.", window, sep = ""))
	## Eliminate all related toolbars
	tkToolDel(paste("$Tk.", window, sep = ""))
	## Possibly delete other resources...
	## ...
	assignTemp(".guiWins", .guiWins)
	tkdestroy(win)  # Actually destroy the window
	invisible(TRUE)
}
