## Menu functions
.menuClear <- function ()
{
	if (.isJGR()) return(invisible(NULL))
	res <- switch(Sys.info()["sysname"],
		Windows = NULL,
		Darwin = .macMenuClear(),
		.unixMenuClear()
	)
	return(invisible(res))
}

.menuFileInit <- function ()
{
		if (.isJGR()) return(invisible(NULL))
		res <- switch(Sys.info()["sysname"],
		Windows = NULL,
		Darwin = NULL, # TODO: should we have a default menu?
		.unixMenuFileInit()
	)
	return(invisible(res))
}

.ctxMenuFileInit <- function ()
{
		res <- switch(Sys.info()["sysname"],
		Windows = NULL,
		Darwin = NULL, # TODO: should we have a default menu?
		.unixCtxMenuFileInit()
	)
	return(invisible(res))
}

.checkMenuName <- function (menuname)
{
	## Make sure menuname is correct...
	menuname <- as.character(menuname)
	if (length(menuname) != 1)
		stop("'menuname' must be a single character string")
	
	## TODO: this is not allowed in JGR!
	## $ConsoleMain/<menu> is equivalent, and thus, transformed into <menu>
	menuname <- sub("^\\$ConsoleMain/", "", menuname)
	
	if (menuname == "")
		stop("You cannot use an empty menuname")
	## Do not accept $ConsoleMain, $ConsolePopup, $Graph<n>Main,
	## or $Graph<n>Popup alone: need a / after the name of a special menu
	if (grepl("^\\$Console(Main|Popup)", menuname) &&
		!grepl("^\\$Console(Main|Popup)/.+", menuname))
		stop("You must define a submenu after the name of a special menu")
	if (grepl("^\\$Graph[0-9]+(Main|Popup)", menuname) &&
		!grepl("^\\$Graph[0-9]+(Main|Popup)/.+", menuname))
		stop("You must define a submenu after the name of a special menu")
	## Return the (possibly arranged) menuname
	return(menuname)
}

menuNames <- function ()
{
	if (.isJGR()) return(.jgrMenuNames())
	res <- switch(Sys.info()["sysname"],
		Windows = winMenuNames(),
		Darwin = .macMenuNames(),
		.unixMenuNames()
	)
	return(res)
}

menuItems <- function (menuname)
{
	menuname <- .checkMenuName(menuname)
	
	if (.isJGR()) return(.jgrMenuItems(menuname))
	res <- switch(Sys.info()["sysname"],
		Windows = winMenuItems(menuname),
		Darwin = .macMenuItems(menuname),
		.unixMenuItems(menuname)
	)
	return(res)
}

menuAdd <- function (menuname)
{
	menuname <- .checkMenuName(menuname)
	
	if (.isJGR()) return(invisible(.jgrMenuAdd(menuname)))
	res <- switch(Sys.info()["sysname"],
		Windows = winMenuAdd(menuname),
		Darwin = .macMenuAdd(menuname),
		.unixMenuAdd(menuname)
	)
	return(invisible(res))
}

menuAddItem <- function (menuname, itemname, action)
{
	menuname <- .checkMenuName(menuname)
	
	if (.isJGR()) return(invisible(.jgrMenuAddItem(menuname, itemname, action)))
	res <- switch(Sys.info()["sysname"],
		Windows = .winMenuAddItem(menuname, itemname, action),
		Darwin = .macMenuAddItem(menuname, itemname, action),
		.unixMenuAddItem(menuname, itemname, action)
	)
	return(invisible(res))
}

menuDel <- function (menuname)
{
	menuname <- .checkMenuName(menuname)
	
	if (.isJGR()) return(invisible(.jgrMenuDel(menuname)))
	res <- switch(Sys.info()["sysname"],
		Windows = try(winMenuDel(menuname), silent = TRUE),
		Darwin = .macMenuDel(menuname),
		.unixMenuDel(menuname)
	)
	return(invisible(res))
}

menuDelItem <- function (menuname, itemname)
{
	menuname <- .checkMenuName(menuname)
	
	if (.isJGR()) return(invisible(.jgrMenuDelItem(menuname, itemname)))
	res <- switch(Sys.info()["sysname"],
		Windows = try(winMenuDelItem(menuname, itemname), silent = TRUE),
		Darwin = .macMenuDelItem(menuname, itemname),
		.unixMenuDelItem(menuname, itemname)
	)
	return(invisible(res))
}


## JGR menus manipulation functions
## The default menus in JGR (for consistence with other implementations, we
## don't want to see them!)
.jgrDefaultMenus <-
	c("File", "Edit", "Workspace", "Packages & Data", "Window", "Help")

.jgrMenuMem <- function ()
{
	## Get an environment with info about JGR menus I cannot get otherwise
	## Basically, I keep track of two things here:
	## 1) which menu action is related to which menu entry
	## 2) for separators, names can be any number of '-', but it is only '-'
	##    in JGR => keep track of the correspondance!
	mnu <- .getTemp(".jgrMenuMem")
	if (is.null(mnu)) {
		mnu <- new.env()
		.assignTemp(".jgrMenuMem", mnu) 
	}
	return(mnu)
}

.jgrMenuMemAdd <- function (menuname, itemname, action)
{
	e <- .jgrMenuMem()
	e[[paste(menuname, itemname, sep = "//")]] <- action
}

.jgrMenuMemGet <- function (menuname, itemname)
{
	e <- .jgrMenuMem()
	e[[paste(menuname, itemname, sep = "//")]]
}

.jgrMenuMemDel <- function (menuname, itemname)
{
	e <- .jgrMenuMem()
	item <- paste(menuname, itemname, sep = "//")
	if (exists(item, envir = e, inherits = FALSE)) rm(list = item, envir = e)	
}

## I redefine the jgr.XXX() function here because I don't want to depend
## on JGR in this package (too much trouble on install of JGR and I never
## want to force users of svDialogs to install JGR)
.jgr.register.function <- function (fun) 
{
    if (is.null(.GlobalEnv$.jgr.user.functions)) 
        .GlobalEnv$.jgr.user.functions <- list()
    fnc <- length(.GlobalEnv$.jgr.user.functions) + 1
    .GlobalEnv$.jgr.user.functions[[fnc]] <- fun
    paste(".jgr.user.functions[[", fnc, "]]()", sep = "")
}

.jgr.getMenuNames <- function () 
{
    if (!.isJGR()) {
        cat(".jgr.getMenuNames() cannot be used outside JGR.\n")
        return(invisible(NULL))
		J <- function (...) return() # Just to avoir R CMD check warning
    }
    J("org/rosuda/JGR/JGR")$getMenuNames()
}

.jgr.getMenuItemNames <- function (menu) 
{
    if (!.isJGR()) {
        cat("/jgr.getMenuItemNames() cannot be used outside JGR.\n")
        return(invisible(NULL))
		J <- function (...) return() # Just to avoir R CMD check warning
    }
    J("org/rosuda/JGR/JGR")$getMenuItemNames(as.character(menu))
}

.jgr.addMenu <- function (name) 
{
    if (!.isJGR()) {
        cat(".jgr.addMenu() cannot be used outside JGR.\n")
        return(invisible(NULL))
		.jcall <- function (...) return() # Just to avoir R CMD check warning
    }
    invisible(.jcall("org/rosuda/JGR/JGR", "V", "addMenu", as.character(name)))
}

.jgr.insertMenu <- function (name, index) 
{
    if (!.isJGR()) {
        cat(".jgr.insertMenu() cannot be used outside JGR.\n")
        return(invisible(NULL))
		.jcall <- function (...) return() # Just to avoir R CMD check warning
    }
    invisible(.jcall("org/rosuda/JGR/JGR", "V", "insertMenu", 
        as.character(name), as.integer(index - 1)))
}

.jgr.addMenuItem <- function (menu, name, command, silent = TRUE) 
{
    if (!.isJGR()) {
        cat(".jgr.addMenuItem() cannot be used outside JGR.\n")
        return(invisible(NULL))
		.jcall <- function (...) return() # Just to avoir R CMD check warning
    }
    if (is.function(command)) 
        command <- .jgr.register.function(command)
    invisible(.jcall("org/rosuda/JGR/JGR", "V", "addMenuItem", 
        as.character(menu), as.character(name), as.character(command), 
        as.logical(silent)))
}

.jgr.insertMenuItem <- function (menu, name, command, index, silent = TRUE) 
{
    if (!.isJGR()) {
        cat(".jgr.insertMenuItem() cannot be used outside JGR.\n")
        return(invisible(NULL))
		.jcall <- function (...) return() # Just to avoir R CMD check warning
    }
    if (is.function(command)) 
        command <- .jgr.register.function(command)
    invisible(.jcall("org/rosuda/JGR/JGR", "V", "insertMenuItem", 
        as.character(menu), as.character(name), as.character(command), 
        as.logical(silent), as.integer(index - 1)))
}

.jgr.addMenuSeparator <- function (menu) 
{
    if (!.isJGR()) {
        cat(".jgr.addMenuSeparator() cannot be used outside JGR.\n")
        return(invisible(NULL))
		.jcall <- function (...) return() # Just to avoir R CMD check warning
    }
    invisible(.jcall("org/rosuda/JGR/JGR", "V", "addMenuSeparator", 
        as.character(menu)))
}

.jgr.insertMenuSeparator <- function (menu, index) 
{
    if (!.isJGR()) {
        cat(".jgr.insertMenuSeparator() cannot be used outside JGR.\n")
        return(invisible(NULL))
		.jcall <- function (...) return() # Just to avoir R CMD check warning
    }
    invisible(.jcall("org/rosuda/JGR/JGR", "V", "insertMenuSeparator", 
        as.character(menu), as.integer(index - 1)))
}

.jgr.addSubMenu <- function (menu, subMenuName, labels, commands) 
{
    if (!.isJGR()) {
        cat(".jgr.addSubMenu() cannot be used outside JGR.\n")
        return(invisible(NULL))
		J <- function (...) return() # Just to avoir R CMD check warning
    }
    invisible(J("org/rosuda/JGR/JGR")$addSubMenu(menu, subMenuName, 
        labels, commands))
}

## There seems to be a bug in the original jgr.insertSubMenu() function
.jgr.insertSubMenu <- function (menu, subMenuName, labels, commands, index) 
{
    if (!.isJGR()) {
        cat(".jgr.addSubMenu() cannot be used outside JGR.\n")
        return(invisible(NULL))
		J <- function (...) return() # Just to avoir R CMD check warning
    }
    invisible(J("org/rosuda/JGR/JGR")$insertSubMenu(menu, subMenuName, 
        as.integer(index - 1), labels, commands))
}

.jgr.removeMenu <- function (index) 
{
    if (!.isJGR()) {
        cat(".jgr.removeMenu() cannot be used outside JGR.\n")
        return(invisible(NULL))
		J <- function (...) return() # Just to avoir R CMD check warning
    }
    J("org/rosuda/JGR/JGR")$removeMenu(as.integer(index - 1))
}

.jgr.removeMenuItem <- function (menu, index) 
{
    if (!.isJGR()) {
        cat(".jgr.removeMenuItem() cannot be used outside JGR.\n")
        return(invisible(NULL))
		J <- function (...) return() # Just to avoir R CMD check warning
    }
    J("org/rosuda/JGR/JGR")$removeMenuItem(as.character(menu), 
        as.integer(index - 1))
}

.jgrRun <- function (cmd, envir = .GlobalEnv)
{
	## This function is used for submenus where it is not currently possible
	## to declare silent = FALSE
	## Basically, it prints the command, then, evaluate it, capturing output
	## and finally, it prints that output
	cat("> ", deparse(substitute(cmd)), "\n", sep = "")
	res <- capture.output(cmd)
	cat(res)
	cat("\n")
}

## Prepare a command to be printed more or less correctly on screen
## for JGR submenus
.jgrAction <- function (cmd) paste("svDialogs:::.jgrRun(", cmd, ")", sep = "")

## Implementation of JGR menus manipulation
.jgrMenuNames <- function ()
{
	## For consistency with the other implementations, do not return
	## the stadard menus File, Edit, Workspace, Packages & Data, Window, Help
	res <- .jgr.getMenuNames()
	res <- res[!res %in% .jgrDefaultMenus]
	return(res)
}

.jgrMenuItems <- function (menuname)
{
	## For consistency with the other implementations, return character(0)
	res <- try(.jgr.getMenuItemNames(menuname), silent = TRUE)
	## if the menu is not found
	if (inherits(res, "try-error")) return(character(0)) else return(res)
}

.jgrMenuAdd <- function (menuname)
{
	## In JGR, one accepts only menus and one level of submenus... and
	## submenus are managed in a quite different way!
	## For compatibility, we allow to use menu/submenu
	mnu <- strsplit(menuname[1], "/", fixed = TRUE)[[1]]
	l <- length(mnu)
	if (l == 0) return(invisible(NULL))
		
	.addTopMenu <- function (topmenu) {
		if (!topmenu %in% .jgr.getMenuNames()) # Here, we check all JGR menus!
		.jgr.addMenu(topmenu)
	}
		
	if (l == 1) {
		.addTopMenu(mnu)
		return(invisible(NULL))
	}
	if (l > 2) stop("Only one submenu level allowed on JGR")
	## Make sure the topmenu is define, and add an empty submenu to it
	.addTopMenu(mnu[1])
	.jgr.addSubMenu(mnu[1], mnu[2], character(0), character(0))
	return(invisible(NULL))
}

.jgrMenuAddItem <- function (menuname, itemname, action)
{
	## In JGR, one accepts only menus and one level of submenus... and
	## submenus are managed in a quite different way!
	## For compatibility, we allow to use menu/submenu
	mnu <- strsplit(menuname[1], "/", fixed = TRUE)[[1]]
	l <- length(mnu)
	if (l == 0) return(invisible(NULL)) # Nothing to do...
		
	silent <- FALSE
	## Special cases for "none", no action associated with this menu
	if (action == "none") {
		action <- ""
		silent <- TRUE
	}
	## Replace \n by \\n, and \t by \\t
	action <- gsub("\n", "\\\\n", action)
	action <- gsub("\t", "\\\\t", action)
		
	if (l == 1) {
		## First, make sure the menu exists
		.jgrMenuAdd(mnu)
		## Are we trying to add a separator?
		if (grepl("^-+$", itemname)) { # This must be a separator
			.jgr.addMenuSeparator(mnu)
		} else { # This must be a menu entry
			## Is the menu entry already implemented?
			items <- .jgrMenuItems(mnu)
			if (itemname %in% items) {
				## Delete and recreate it with the new action
				idx <- (1:length(items))[items == itemname][1]
				.jgr.removeMenuItem(mnu, idx)
				if (action == "enable")
					action <- .jgrMenuMemGet(mnu, itemname)
				if (is.null(action)) action <- ""
				if (action == "disable") {
					action <- 'cat("- disabled menu item...\n")'
					silent <- TRUE
					.jgr.insertMenuItem(mnu, itemname, action, idx,
						silent = silent)
					## Don't change action in .jgrMenus() so that we can
					## recover it with "enable"!
				} else {
					.jgr.insertMenuItem(mnu, itemname, action, idx,
						silent = silent)	
					.jgrMenuMemAdd(mnu, itemname, action)
				}
			} else {
				## This is a new item => just add it
				if (action == "enable") action <- ""
				if (action == "disable") {
					action <- 'cat("- disabled menu item...\n")'
					silent <- TRUE
					.jgr.addMenuItem(mnu, itemname, action, silent = silent)
					.jgrMenuMemAdd(mnu, itemname, "")
				} else {
					.jgr.addMenuItem(mnu, itemname, action, silent = silent)
					.jgrMenuMemAdd(mnu, itemname, action)
				}
			}
		}			
		return(invisible(NULL))
	}
		
	if (l == 2) {
		## We add an entry in a submenu. In JGR, we must delete and
		## reconstruct the submenu entirely!
		## Does this submenu already exists?
		## (note: in JGR it can be a menu entry as well!)
		items <- .jgrMenuItems(mnu[1])
		if (mnu[2] %in% items) { # The submenu already exists...
			## Delete it and reconstruct it with the added or changed item
			idx <- (1:length(items))[items == mnu[2]][1]
			.jgr.removeMenuItem(mnu[1], idx)
			## Get the list of entries in the submenu from .jgrMenus
			## (how to get it otherwise???)
			actions <- .jgrMenuMemGet(mnu[1], mnu[2])
			if (is.null(actions)) { # Apparently nothing in there yet
				if (action == "enable") action <- ""
				if (action == "disable") {
					action <- 'cat("- disabled menu item...\n")'
					.jgr.insertSubMenu(mnu[1], mnu[2], c(itemname, "-"),
						c(action, ""), idx)
					actions <- ""
				} else {
					action <- .jgrAction(action)
					.jgr.insertSubMenu(mnu[1], mnu[2], c(itemname, "-"),
						c(action, ""), idx)
					actions <- action
				}
				names(actions) <- itemname
				.jgrMenuMemAdd(mnu[1], mnu[2], actions)
			} else { # There are already items in this submenu
				if (action == "disable") { # We want to disable one action
					if (!itemname %in% names(actions))
						return(invisible(NULL))
					actions[[itemname]] <- 'cat("- disabled menu item...\n")'
					if (length(actions) == 1) {
						.jgr.insertSubMenu(mnu[1], mnu[2],
							c(names(actions), "-"), c(actions, ""), idx)	
					} else {
						.jgr.insertSubMenu(mnu[1], mnu[2],
							names(actions), actions, idx)
					}
				} else if (action == "enable") {
					if (!itemname %in% names(actions))
						return(invisible(NULL))
					if (length(actions) == 1) {
						.jgr.insertSubMenu(mnu[1], mnu[2],
							c(names(actions), "-"), c(actions, ""), idx)
					} else {
						.jgr.insertSubMenu(mnu[1], mnu[2],
							names(actions), actions, idx)	
					}
				} else { # An action is defined
					actions[[itemname]] <- .jgrAction(action)
					if (length(actions) == 1) {
						.jgr.insertSubMenu(mnu[1], mnu[2],
							c(names(actions), "-"), c(actions, ""), idx)
					} else {
						.jgr.insertSubMenu(mnu[1], mnu[2],
							names(actions), actions, idx)	
					}
					.jgrMenuMemAdd(mnu[1], mnu[2], actions)	
				}
			}
		} else { # The submenu does not exists yet, create it now
			## First, make sure the top menu exists
			.jgrMenuAdd(mnu[1])
			if (action == "enable") action <- ""
			if (action == "disable") {
				action <- 'cat("- disabled menu item...\n")'
				.jgr.addSubMenu(mnu[1], mnu[2], c(itemname, "-"), c(action, ""))
				actions <- ""
			} else {
				action <- .jgrAction(action)
				.jgr.addSubMenu(mnu[1], mnu[2], c(itemname, "-"), c(action, ""))
				actions <- action
			}
			names(actions) <- itemname
			.jgrMenuMemAdd(mnu[1], mnu[2], actions)
		}
		return(invisible(NULL))
	}
	
	if (l > 2) stop("Only one submenu level allowed on JGR")
}

.jgrMenuDel <- function (menuname)
{
	## In JGR, one accepts only menus and one level of submenus... and
	## submenus are managed in a quite different way!
	## For compatibility, we allow to use menu/submenu
	mnu <- strsplit(menuname[1], "/", fixed = TRUE)[[1]]
	l <- length(mnu)
	if (l == 0) return(invisible(NULL))
		
	if (l == 1) {
		if (mnu %in% .jgrDefaultMenus) # Do not allow to delete default menus
			return(invisible(NULL))
		## Get the position of this menu and make sure it is not a default menu!
		allmnu <- .jgr.getMenuNames()
		allpos <- 1:length(allmnu)
		pos <- allpos[allmnu == mnu]
		if (!length(pos)) return(invisible(NULL)) # Not found
		pos <- rev(pos)[1]
		.jgr.removeMenu(pos)
	}
	if (l == 2) # Remove a submenu (note, that for JGR, this is the same as removing a menu item!)
		.jgrMenuDelItem(mnu[1], mnu[2])
	## If there are more levels, do nothing because these submenus do not exist!
	return(invisible(NULL))
}

.jgrMenuDelItem <- function (menuname, itemname)
{
	## On JGR, all separators are named "-", but on, e.g., Windows, I must
	## use a different name for each separator => What should we do???
	## Here, we consider '-' for the first one, '--' for the second one, etc.
	mnu <- strsplit(menuname[1], "/", fixed = TRUE)[[1]]
	l <- length(mnu)
	if (l == 0) return(invisible(NULL))
		
	if (l == 1) {
		## Are we trying to delete a separator?
		items <- .jgrMenuItems(mnu)
		if (grepl("^-+$", itemname)) { # This must be a separator
			isSep <- items == "-"
			posSep <- (1:length(items))[isSep]
			## Depending on the number of minus signs we try to remove
			## first, second, etc. separator
			nsep <- nchar(itemname)
			if (length(posSep) < nsep) return(invisible(NULL))
			idx <- posSep[nsep]
		} else { # This must be a menu entry
			idx <- (1:length(items))[itemname == items]
		}
		if (!length(idx)) return(invisible(NULL))
		.jgr.removeMenuItem(mnu, idx)
		.jgrMenuMemDel(mnu, itemname)
	}
		
	if (l == 2) {
		## We want to eliminate an item from a submenu. In JGR, there is no
		## function for that, but we can delete and recreate the submenu
		## without this item!
		items <- .jgrMenuItems(mnu[1])
		if (!mnu[2] %in% items) return(invisible(NULL)) # Submenu not there
		## Get the list of submenus currently defined from our cache version
		actions <- .jgrMenuMemGet(mnu[1], mnu[2])
		if (is.null(actions) || !itemname %in% names(actions)) # Apparently not there
			return(invisible(NULL))
		## Delete and reconstruct the submenu without itemname
		idx <- (1:length(items))[items == mnu[2]][1]
		.jgr.removeMenuItem(mnu[1], idx)
		## Recreate the submenu after eliminating itemname
		actions <- actions[names(actions) != itemname]
		if (!length(actions)) {
			.jgr.insertSubMenu(mnu[1], mnu[2], character(0), character(0), idx)
		} else if (length(actions) == 1) {
			.jgr.insertSubMenu(mnu[1], mnu[2], c(names(actions), "-"),
				c(actions, ""), idx)
		} else {
			.jgr.insertSubMenu(mnu[1], mnu[2], names(actions), actions, idx)
		}
		.jgrMenuMemAdd(mnu[1], mnu[2], actions)
	}
	return(invisible(NULL)) # Not more than one submenu level allowed!  
}


## Windows version and standard winMenuXXX
## TODO: fallback system for Rterm???
.winMenuAddItem <- function (menuname, itemname, action)
{
	## As in R 2.14.1, the original winMenuAddItem() does things I don't like
	## much when using 'enable' or 'disable' for action, on a non existing menu:
	## it creates it with the action being 'enable' or 'disable'... I suppose
	## if is a feature, but I want a different behaviour here: to ignore such
	## a command applied to a non-existing menu item!
	if (action %in% c("enable", "disable")) {
		menus <- winMenuItems(menuname)
		if (!is.null(menus) && itemname %in% names(menus)) {
			## The menu exists... enable or disable it!
			return(winMenuAddItem(menuname, itemname, action))
		} else return(invisible(NULL))
	} else return(winMenuAddItem(menuname, itemname, action))
}

## Mac OS X version
## Note, either we use AppleScript folder (by default) or the XMenu folder
## Install XMenu from http://xmenu.en.softonic.com/mac.
## You should use the custom commands only for R, because it will erase
## everything in it everytime the svDialogs package starts!
## Configure XMenu to display only User-Defined items, and name it R.
## Select folders before files, for icons: None, Big font
## and for menu titles: Text and then, Start at login
#options(useXMenu = TRUE)
.macMenuFolder <- function ()
{
	## Get the root folder for the R menus, depends on wether we use XMenu or not
	useXMenu <- getOption("useXMenu", default = NULL)
	if (is.null(useXMenu)) {
		## If not specified, look if a "R" file or folder exists
		useXMenu <- file.exists("~/Library/Application Support/XMenu/R")
	} else useXMenu <- isTRUE(useXMenu)	
	return(getOption("menuFolder", default = if (useXMenu)
		"~/Library/Application Support/XMenu/Custom" else
		"~/Library/Scripts/Applications/R"))
}

.macMenuClear <- function () {
	#stop("Not implemented yet!")
	
#	## To be called when svDialogs package loads: make sure to zap all
#    ## custom menu items that may have been previously defined
#    ## (also call it when the package closes)
#    odir <- getwd()
#    on.exit(setwd(odir))
#    setwd(.macMenuFolder())
#	setwd("..")
#    folder <- file.path(".", basename(.macMenuFolder()))
#	unlink(folder, recursive = TRUE)
#    dir.create(folder, recursive = TRUE)
#    ## Now, I can assume that the dir is created and is empty
#    return(invisible(NULL))
}

.macMenuNames <- function ()
{
	#stop("Not implemented yet!")
}

.macMenuItems <- function (menuname)
{
	#stop("Not implemented yet!")
}

.macMenuAdd <- function (menuname)
{
    #stop("Not implemented yet!")
	
#	## Menus are folders created in ~/Scripts/Applications/R/Custom
#    ## I just need to create (recursively) the directories
#    dir.create(file.path(.macMenuFolder(), menuname),
#		showWarnings = FALSE, recursive = TRUE)
#    return(invisible(NULL))
}

.macMenuAddItem <- function (menuname, itemname, action)
{
	#stop("Not implemented yet!")
	
#	## TODO: manage 'enable' and 'disable'!!!
#	## Make sure that the dir is created
#    .macMenuAdd(menuname)
#    ## Switch to this folder
#	odir <- getwd()
#    on.exit(setwd(odir))
#    setwd(file.path(.macMenuFolder(), menuname))
#	## Add an executable file in it with 'itemname' name
#	## that contains AppleScript code to run action in R
#	## Determine if R is run in R.app or in a terminal window
#	if (.Platform$GUI == "AQUA") {
#		## Can be R or R64 or SciViews R or SciViews R64!
#		app <- paste('"', system("osascript -e 'name of application \"R\"'",
#			intern = TRUE), '"', sep = "")
#	} else app <- "\"Terminal\""
#	## Define action accordingly
#	if (action == "none") {
#		cmd <- "to activate"
#	} else {
#		## Make sure to quote "
#		action <- gsub('"', '\\\\"', action)
#		## Also replace \n, \r and \t
#		action <- gsub('\n', '\\\\\\\\n', action)
#		action <- gsub('\r', '\\\\\\\\r', action)
#		action <- gsub('\t', '\\\\\\\\t', action)
#		if (app == "\"Terminal\"") {
#			cmd <- paste("to do script \"", action, "\" in window 1", sep = "")	
#		} else {
#			cmd <- paste("to cmd \"", action, "\"", sep = "")
#		}
#	}
#	## Compile applescript item
#	system(paste("osacompile -e 'tell application ", app, " ", cmd,
#		"' -o \"", itemname, ".app\"", sep = ""), ignore.stdout = TRUE,
#		ignore.stderr = TRUE)
#    return(invisible(NULL))
}

.macMenuDel <- function (menuname)
{
	#stop("Not implemented yet!")
	
#	## Unlink does not like ~ => change working dir first
#    odir <- getwd()
#    on.exit(setwd(odir))
#    setwd(.macMenuFolder())
#    unlink(menuname, recursive = TRUE)
#    return(invisible(NULL))
}

.macMenuDelItem <- function (menuname, itemname)
{
	#stop("Not implemented yet!")

#    ## Unlink does not like ~ => change working dir first
#    odir <- getwd()
#    on.exit(setwd(odir))
#    setwd(file.path(.macMenuFolder()))
#	unlink(file.path(".", menuname, paste(itemname, "app", sep = ".")),
#		recursive = TRUE) 
#    return(invisible(NULL))    
}


## This holds the custom menu structure in an R object
.Rmenu <- function ()
{
	## The custom R menu is cached in a Rmenu object in SciViews:TempEnv
	return(.getTemp("Rmenu", default = list(), mode = "list"))
}

## Linux/Unix version
## To use R custom context menu, you have to install xvkbd, xdotool,
## zenity (and, possibly, yad)
## You need also to compile and install ctxmenu
## On Ubuntu:
## sudo apt-get install xvkbd xdotool zenity
## Warning, you need to insstall the English (US) keymap, even if you don't use
## it. Otherwise, xvkbd will issue strange things in your R console!
## TODO: install and configure ctxmenu... + add shortcut keys!
## Use xbindkeys to bind shell commands to keyboard and mouse keys
## chmod +x ctxmenu

## We need to write menu files in /tmp... but CRAN policies do not allow this
## unless we got acknowledgement by the user in the interactive session.
## This function checks that R runs in interactive mode and the user gave
## acknowledgement, either interactively,  or by defining the option
## 'svDialogs.tmpfiles' to TRUE
.tmpfilesAllowed <- function ()
	return(interactive() && isTRUE(getOption("svDialogs.tmpfiles", FALSE)))

.unixTmpfilesAsk <- function ()
{
	if (!interactive()) return(FALSE)
	
	## Make sure the user gave explicit right to create custom menus temp
	## files, either interactively, or through the setting of an option
	## options(svDialogs.tmpfiles = TRUE)
	opt <- getOption("svDialogs.tmpfiles")
	if (is.null(opt)) { # Ask user interactively
		if (okCancelBox("Install custom menu configuration files in ~/.ctxmenu/tmp/?")) {
			options(svDialogs.tmpfiles = TRUE)
			## Make sure to clear old menus, and to install new ones
			.menuClear()
			.menuFileInit()
			.ctxMenuFileInit()	
			## Make sure that ctxmenu is installed
			if (Sys.which("ctxmenu") == "") {
				warning("Menus will not be displayed if you do not install ctxmenu properly, see: http://www.sciviews.org/SciViews-R/ctxmenu.zip")
				return(FALSE)
			} else return(TRUE) # Everything should be ok!
		} else {
			options(svDialogs.tmpfiles = FALSE)
			## Indicate that it will not be possible to use custom menus before
			## allowed by the corresponding option
			warning("Menus will not be displayed unless you agree to create config files using options(svDialogs.tmpfiles = TRUE)")
			return(FALSE)
		}
	} else if (!isTRUE(opt)) {
		warning("Menus will not be displayed unless you agree to create config files using options(svDialogs.tmpfiles = TRUE)")
		return(FALSE)
	} else return(TRUE) # Note: we do not check again for ctxmenu here!
}

.unixMenuFolder <- function ()
{
	## Get the root folder for the R menus
	mnudir <- getOption("menuFolder", default = "~/.ctxmenu/tmp")
	## Make sure this directory exists, in one can write to it!
	if (.tmpfilesAllowed() && !file.exists(mnudir))
		try(dir.create(mnudir, showWarnings = FALSE, recursive = TRUE),
			silent = TRUE)
	return(mnudir)
}

.unixMenuFile <- function ()
{
	## Get the name of the file that contains the R menu
	winid <- .getTemp(".winid", default = Sys.getenv("WINDOWID"))
	.assignTemp(".winid", winid)
	## Do not use user name in the filename (B. Ripley's request)
	#user <- .getTemp(".user", default = Sys.getenv("USER"))
	#.assignTemp(".user", user)
	#return(file.path(.unixMenuFolder(), paste(user, winid,
	#	"Menu.txt", sep = "")))
	return(file.path(.unixMenuFolder(), paste(winid,
		"Menu.txt", sep = "")))
}

.unixMenuFileInit <- function ()
{
	## Can we generate files in /tmp?
	if (!.tmpfilesAllowed()) return(invisible(NULL))
	
	## Initialize the R menu file with default items
	fil <- .unixMenuFile()
	## Get the default R menu and start from there
	def <- getOption("RMenuFile",
		default = file.path("~", ".ctxmenu", "RMenu.txt"))
	if (file.exists(def)) {
		file.copy(def, fil, overwrite = TRUE)
	} else file.copy(system.file("gui", "RMenuLinux.txt",
		package = "svDialogs"), fil, overwrite = TRUE)
	return(invisible(NULL))
}

.unixCtxMenuFile <- function ()
{
	## Get the name of the file that contains the R context menu
	winid <- .getTemp(".winid", default = Sys.getenv("WINDOWID"))
	.assignTemp(".winid", winid)
	## Do not use user name in the filename (B. Ripley's request)
	#user <- .getTemp(".user", default = Sys.getenv("USER"))
	#.assignTemp(".user", user)
	#return(file.path(.unixMenuFolder(), paste(user, winid,
	#	"CtxMenu.txt", sep = "")))
	return(file.path(.unixMenuFolder(), paste(winid,
		"CtxMenu.txt", sep = "")))
}

.unixCtxMenuFileInit <- function ()
{
	## Can we generate files in /tmp?
	if (!.tmpfilesAllowed()) return(invisible(NULL))
	
	## Initialize the R context menu file with default items
	fil <- .unixCtxMenuFile()
	## Get the default R context menu and start from there
	def <- getOption("RCtxMenuFile",
		default = file.path("~", ".ctxmenu", "RCtxMenu.txt"))
	if (file.exists(def)) {
		file.copy(def, fil, overwrite = TRUE)
	} else file.copy(system.file("gui", "RCtxMenuLinux.txt",
		package = "svDialogs"), fil, overwrite = TRUE)
	return(invisible(NULL))
}

.unixMenuSave <- function (mnu, file = TRUE)
{
	## Save the menu structure in Rmenu object in SciViews:TempEnv and in a file
	## mnu is either a list of lists with menu entries, or NULL to delete all
	## custom menus
	.assignTemp("Rmenu", mnu)
	
	# Do nothing on files, unless interactive() and got user's acceptation
	if (!.tmpfilesAllowed()) return(invisible(NULL))
	if (!isTRUE(file)) return(invisible(NULL))
	## The menu file is:
	fil <- .unixMenuFile()
	ctxfil <- .unixCtxMenuFile()
	if (is.null(mnu)) {
		## Clear the file
		unlink(fil)
	} else {
		## Populate the file with the content of the Rmenu object
		makeMenu <- function (lst, indent = 0, file = fil, ctxfile = ctxfil) {
			l <- length(lst)
			if (l < 1) return()
			nms <- names(lst)
			for (i in 1:l) {
				item <- nms[i]
				if (is.list(lst[[i]])) {
					## Special case for '$ConsolePopup'
					if (item == "$ConsolePopup" && !is.null(ctxfile)) {
						makeMenu(lst[[i]], indent = 0, file = ctxfile, ctxfile = NULL)
					} else {
						## Create a new menu
						cat("\n", rep("\t", indent), "submenu=", item, "\n",
							sep = "", file = file, append = TRUE)
						makeMenu(lst[[i]], indent = indent + 1, file = file, ctxfile = NULL)
					}
				} else {
					## Is this a separator?
					if (grepl("^-",  item)) {
						cat("\n", rep("\t", indent), "separator\n",
							sep = "", file = file, append = TRUE)
					} else { # Add an item in current menu
						ind <- rep("\t", indent)
						## Rework commands using xvkbd -text "cmd\r"
						cmd <- as.character(lst[[i]])[1]
						if (cmd == "none" || !is.null(attr(lst[[i]], "state"))) {
							cmd <- "NULL" # This is the "no cmd" or "disabled" for ctxmenu
						} else {
							cmd <- paste(cmd, "\\n", sep = "")
							cmd <- paste("xvkbd -text", shQuote(cmd))
						}
						cat("\n", ind, "item=", item, "\n", ind, "cmd=", cmd,
							"\n", sep = "", file = file, append = TRUE)
					}
				}
			}
		}
		## Initialize the menu with default items
		.unixMenuFileInit()
		.unixCtxMenuFileInit()
		## Add custom menus to it...
		cat("\nseparator # Here starts the custom R menus\n\n",
			file = fil, append = TRUE)
		cat("\nseparator # Here starts the custom R context menus\n\n",
			file = ctxfil, append = TRUE)
		makeMenu(mnu)
	}
	return(invisible(fil))
}

.unixMenuClear <- function ()
{
    ## To be called when svDialogs package loads: make sure to zap all
    ## custom menu items that may have been previously defined
    ## (also call it when the package closes)

	# Do nothing on files, unless interactive() and got user's acceptation
	if (!.tmpfilesAllowed()) return(invisible(NULL))
	
	## Also clear the local object
	.unixMenuSave(NULL)
	unlink(.unixMenuFile())
	unlink(.unixCtxMenuFile())
    return(invisible(NULL))
}

.unixMenuNames <- function (mnu = .Rmenu(), parent = character(0))
{
	if (!length(mnu)) return(character(0))
	## List all custom menu and (sub)menu
	## Iteratively traverse the list and enumerate all sublists (= submenus)
	items <- names(mnu)
	submnus <- character(0)
	for (item in items)
		if (is.list(mnu[[item]]))
			submnus <- c(submnus, paste(parent, item, sep = ""),
				.unixMenuNames(mnu[[item]],
					parent = paste(parent, item, "/", sep = "")))
	names(submnus) <- NULL
	## Eliminate $ConsolePopup, $Graph<n>Main & $Graph<n>Popup
	submnus <- submnus[submnus != "$ConsolePopup"]
	isGraphMenu <- grepl("^\\$Graph[0-9]+(Main|Popup)$", submnus)
	submnus <- submnus[!isGraphMenu]
	return(submnus)
}

.unixMenuItems <- function (menuname)
{
	## List all menu items in a given (sub)menu
	mnu <- .Rmenu()
	items <- strsplit(as.character(menuname), "/", fixed = TRUE)[[1]]
	## Traverse the menu hierarchy
	for (i in items) {
		mnu <- mnu[[i]]
		if (is.null(mnu) || !is.list(mnu))
			stop("unable to retrieve items for ", menuname,
				" (menu does not exist)")
	}
	if (length(mnu) == 0) return(character(0))
	## Set all submenu items to NULL
	for (i in 1:length(mnu))
		if (is.list(mnu[[i]])) mnu[[i]] <- NULL
	if (length(mnu) == 0) return(character(0)) else return(unlist(mnu))
}

.unixMenuAdd <- function (menuname, itemname = NULL, action = "none") {
	## Make sure we can install required files in ~/.ctxmenu/tmp/
	.unixTmpfilesAsk()
	
	## Add this menu to our Rmenu object
	mnu <- .Rmenu()
	items <- strsplit(as.character(menuname), "/", fixed = TRUE)[[1]]
	## Allow for a maximum of 5 sublevels (should be largely enough!)
	l <- length(items)
	if (l == 1) {
		if (!is.null(mnu[[items[1]]])) {
			## If this is not a list, we got an error
			if (!is.list(mnu[[items[1]]]))
				stop(menuname, " is already defined and is not a menu")
		} else { # Create it
			mnu[[items[1]]] <- list()
		}
		## Do we create an menu item there too?
		if (!is.null(itemname))
			mnu[[items[1]]][[itemname]] <- action
	} else if (l == 2) {
		if (!is.null(mnu[[items[1]]][[items[2]]])) {
			## If this is not a list, we got an error
			if (!is.list(mnu[[items[1]]][[items[2]]]))
				stop(menuname, " is already defined and is not a menu")
		} else { # Create it
			mnu[[items[1]]][[items[2]]] <- list()
		}
		## Do we create an menu item there too?
		if (!is.null(itemname))
			mnu[[items[1]]][[items[2]]][[itemname]] <- action
	} else if (l == 3) {
				if (!is.null(mnu[[items[1]]][[items[2]]][[items[3]]])) {
			## If this is not a list, we got an error
			if (!is.list(mnu[[items[1]]][[items[2]]][[items[3]]]))
				stop(menuname, " is already defined and is not a menu")
		} else { # Create it
			mnu[[items[1]]][[items[2]]][[items[3]]] <- list()
		}
		## Do we create an menu item there too?
		if (!is.null(itemname))
			mnu[[items[1]]][[items[2]]][[items[3]]][[itemname]] <- action
	} else if (l == 4) {
		if (!is.null(mnu[[items[1]]][[items[2]]][[items[3]]][[items[4]]])) {
			## If this is not a list, we got an error
			if (!is.list(mnu[[items[1]]][[items[2]]][[items[3]]][[items[4]]]))
				stop(menuname, " is already defined and is not a menu")
		} else { # Create it
			mnu[[items[1]]][[items[2]]][[items[3]]][[items[4]]] <- list()
		}
		## Do we create an menu item there too?
		if (!is.null(itemname))
			mnu[[items[1]]][[items[2]]][[items[3]]][[items[4]]][[itemname]] <- action
	} else if (l == 5) {
		if (!is.null(mnu[[items[1]]][[items[2]]][[items[3]]][[items[4]]][[items[5]]])) {
			## If this is not a list, we got an error
			if (!is.list(mnu[[items[1]]][[items[2]]][[items[3]]][[items[4]]][[items[5]]]))
				stop(menuname, " is already defined and is not a menu")
		} else { # Create it
			mnu[[items[1]]][[items[2]]][[items[3]]][[items[4]]][[items[5]]] <- list()
		}
		## Do we create an menu item there too?
		if (!is.null(itemname))
			mnu[[items[1]]][[items[2]]][[items[3]]][[items[4]]][[items[5]]][[itemname]] <- action
	} else if (l > 5) {
		stop("You cannot use more than 5 menu levels")
	}
	## Save these changes
	.unixMenuSave(mnu)
	return(invisible(NULL))
}

.unixMenuAddItem <- function (menuname, itemname, action) {
	## Make sure we can install required files in ~/.ctxmenu/tmp/
	.unixTmpfilesAsk()
	
	if (action %in% c("enable", "disable")) {
		## Enable or disable an existing menu item
		if (action == "enable") action <- NULL # To eliminate the attribute
		mnu <- .Rmenu()
		items <- strsplit(as.character(menuname), "/", fixed = TRUE)[[1]]
		## allow for a maximum of 5 sublevels (should be largely enough!)
		l <- length(items)
		if (l == 1) {
			if (!is.null(mnu[[items[1]]][[itemname]]))
				attr(mnu[[items[1]]][[itemname]], "state") <- action
		} else if (l == 2) {
			if (!is.null(mnu[[items[1]]][[items[2]]][[itemname]]))
				attr(mnu[[items[1]]][[items[2]]][[itemname]], "state") <- action
		} else if (l == 3) {
			if (!is.null(mnu[[items[1]]][[items[2]]][[items[3]]][[itemname]]))
				attr(mnu[[items[1]]][[items[2]]][[items[3]]][[itemname]], "state") <- action
		} else if (l == 4) {
			if (!is.null(mnu[[items[1]]][[items[2]]][[items[3]]][[items[4]]][[itemname]]))
				attr(mnu[[items[1]]][[items[2]]][[items[3]]][[items[4]]][[itemname]], "state") <- action
		} else if (l == 5) {
			if (!is.null(mnu[[items[1]]][[items[2]]][[items[3]]][[items[4]]][[items[5]]][[itemname]]))
				attr(mnu[[items[1]]][[items[2]]][[items[3]]][[items[4]]][[items[5]]][[itemname]], "state") <- action
		} else if (l > 5) {
			stop("You cannot use more than 5 menu levels")
		}
		## Save these changes
		.unixMenuSave(mnu)
		return(invisible(NULL))
	} else return(.unixMenuAdd(menuname, itemname, action))
}

.unixMenuDel <- function (menuname) {
	mnu <- .Rmenu()
	items <- strsplit(as.character(menuname), "/", fixed = TRUE)[[1]]
	## Allow for a maximum of 5 sublevels (should be largely enough!)
	l <- length(items)
	if (l == 1 && !is.null(mnu[[items[1]]])) {
		mnu[[items[1]]] <- NULL
	} else if (l == 2 && !is.null(mnu[[items[1]]][[items[2]]])) {
		mnu[[items[1]]][[items[2]]] <- NULL
	} else if (l == 3 && !is.null(mnu[[items[1]]][[items[2]]][[items[3]]])) {
		mnu[[items[1]]][[items[2]]][[items[3]]] <- NULL
	} else if (l == 4 && !is.null(mnu[[items[1]]][[items[2]]][[items[3]]][[items[4]]])) {
		mnu[[items[1]]][[items[2]]][[items[3]]][[items[4]]] <- NULL
	} else if (l == 5 && !is.null(mnu[[items[1]]][[items[2]]][[items[3]]][[items[4]]][[items[5]]])) {
		mnu[[items[1]]][[items[2]]][[items[3]]][[items[4]]][[items[5]]] <- NULL
	} else return(invisible(NULL))
	## Save these changes
	.unixMenuSave(mnu)
	return(invisible(NULL))	
}

.unixMenuDelItem <- function (menuname, itemname) {
	mnu <- .Rmenu()
	items <- strsplit(as.character(menuname), "/", fixed = TRUE)[[1]]
	## Allow for a maximum of 5 sublevels (should be largely enough!)
	l <- length(items)
	if (l == 1 && !is.null(mnu[[items[1]]][[itemname]])) {
		mnu[[items[1]]][[itemname]] <- NULL
	} else if (l == 2 && !is.null(mnu[[items[1]]][[items[2]]][[itemname]])) {
		mnu[[items[1]]][[items[2]]][[itemname]] <- NULL
	} else if (l == 3 && !is.null(mnu[[items[1]]][[items[2]]][[items[3]]][[itemname]])) {
		mnu[[items[1]]][[items[2]]][[items[3]]][[itemname]] <- NULL
	} else if (l == 4 && !is.null(mnu[[items[1]]][[items[2]]][[items[3]]][[items[4]]][[itemname]])) {
		mnu[[items[1]]][[items[2]]][[items[3]]][[items[4]]][[itemname]] <- NULL
	} else if (l == 5 && !is.null(mnu[[items[1]]][[items[2]]][[items[3]]][[items[4]]][[items[5]]][[itemname]])) {
		mnu[[items[1]]][[items[2]]][[items[3]]][[items[4]]][[items[5]]][[itemname]] <- NULL
	} else return(invisible(NULL))
	## Save these changes
	.unixMenuSave(mnu)
	return(invisible(NULL))	  
}
