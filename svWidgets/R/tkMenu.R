tkMenuItems <- function (menu)
{
	M <- getTemp(".guiMenus")[[menu]]
	if (is.null(M))
		stop("unable to retrieve items for ", menu, "\n(menu ", menu,
			" does not exist)")
	Mitems <- M$Items
	if (is.null(Mitems)) {
		res <- character(0)
	} else {
		res <- as.character(Mitems$action)
	}
	if (length(res) > 0) names(res) <- Mitems$name
	res
}

tkMenuAdd <- function (menu, tearoff = FALSE)
{
	## Get the menu name
	Mname <- basename(menu)
	## Get the name of the parent
	Pname <- dirname(menu)
	## Look if the parent exists (must be in .guiMenus in SciViews:TempEnv)
	.guiMenus <- getTemp(".guiMenus")
	## If .guiMenus was not there, create it
	if (is.null(.guiMenus)) {
		.guiMenus <- list()
		class(.guiMenus) <- c("guiMenu", "gui", class(.guiMenus))
	}
	## Do not create the menu if it already exists
	if (menu %in% names(.guiMenus))
		stop("menu ", menu, " already exists!")
	Parent <- .guiMenus[[Pname]]
	if (is.null(Parent)) {
        ## If base menu is a "root", try to create the corresponding Tk top menu
        if (regexpr("/", Pname) < 0) {
			## Get the name of the Tk window that hosts the menu
			TkWinName <- sub("^[$]Tk[.]", "", Pname)
			## Look if such a Tk window is recorded in .gui.Wins
			TkWin <- winGet(TkWinName)
			if (is.null(TkWin))
				stop("menu does not exist, ",
					 "and the parent Tk window cannot be found!")
			## Create the menu
			Parent <- tkmenu(TkWin)
			tkconfigure(TkWin, menu = Parent)
			## Add an entry for this top menu in .guiMenus
			.guiMenus[[Pname]] <- Parent
		} else stop("unable to add menu\n(base menu does not exists)")
	}
	## Check that the name is not already used (for a menu entry, for instance)
    items <- Parent$Items
    if (!is.null(items) && Mname %in% items$name)
		stop("the menu name is already used in this menu!")
	## Now that the parent menu exists, we can create this menu
	## Look where to place the underline (menu shortcut)
    Under <- regexpr("&", Mname)
	if (Under < 0) {  # No '&', place the underline at first position
		Under <- 0
		item <- Mname
	} else {
		Under <- Under - 1  # Because Tk starts numbering at 0
		item <- sub("&", "", Mname)
	}
	Child <- tkmenu(Parent, tearoff = tearoff)
    tkadd(Parent, "cascade", label = item, menu = Child, underline = Under)
 	## Add an entry for this child menu in .guiMenus
	.guiMenus[[menu]] <- Child
    ## ... and register it in the items of parent menu in .guiMenus
	if (tearoff) options <- "tearoff = TRUE" else options <- "tearoff = FALSE"
	entry <- data.frame(name = I(Mname), action = I("[menu]"), image = (""),
		accel = I(""), options = I(options))
	if (is.null(items)) {
		.guiMenus[[Pname]]$Items <- entry
	} else {
        .guiMenus[[Pname]]$Items <- rbind(items, entry)
	}
	## Update the SciViews:TempEnv version of .guiMenus
	assignTemp(".guiMenus", .guiMenus)
	invisible(menu)
}

tkMenuItemCall <- function (expr) {
	## Remove { and } from the deparsed expression
	text <- head(deparse(substitute(expr))[-1], -1)
	cat(">> ", text, "\n")
	timestamp(text, prefix = "", suffix = "") # Was: .Internal(addhistory(text))
	eval(expr, envir = parent.frame())
}

tkMenuAddItem <- function (menu, item, action, image = "", accel = "",
options = "")
{	
	## Look if the menu exists (must be in .guiMenus in SciViews:TempEnv)
	.guiMenus <- getTemp(".guiMenus")
	M <- .guiMenus[[menu]]
	if (is.null(M)) {
		## On the contrary to winMenuAddItem(), if the menu does not exist yet
		## generate an error (but we can later change this behaviour, of course!)
		stop("menu does not exist!")
	}
	## Look if the menu already exists
	Items <- M$Items
	if (!is.null(Items) && item %in% Items$name) {
        ## On the contrary to winMenuAddItem(), it is not allowed to add twice
		## the same menu (would change value, but in Tk, it would add a second
		## time the same menu)
		stop("the menu item ", item, " already exists!")
	}
	## Add the entry at the end of the Tk menu (### TODO: allow other positions)
	## First look if it is a command or a separator
	if (regexpr("^-+$", item) > 0) { #This must be a command
		tkadd(M, "separator")
		action <- "[separator]"
		options <- ""
	} else {  # This is a menu command
		## Look for the '&', indicating where the underline should be located
		Under <- regexpr("&", item)
		if (Under < 0) {  # No '&', place the underline at first position
			Uopt <- ", underline = 0"
			lbl <- item
		} else {
			Uopt <- paste(", underline =", Under - 1) # Tk starts numbering at 0
			lbl <- sub("&", "", item)
		}
		## Do we have to add an image to the menu?
		if (image != "") {
			## Look if the image resource is available
			Img <- imgGet(image)
			if (!is.null(Img)) {
				Iopt <- paste(", image = '", as.character(Img), "'", sep ="")
			} else Iopt <- ""
		} else Iopt <- ""
        ## Do we have an accelerator defined for this menu?
		if (accel != "") {
			Aopt <- paste(', accelerator = "', accel, '"', sep = "")
			## Compute the Tk accelerator and make corresponding binding
			tkAccel <- paste("<", tolower(accel), ">", sep ="")
			## 'ctrl+' becomes 'Control-'
			tkAccel <- sub("ctrl[+]", "Control-", tkAccel)
			## 'shift+' becomes 'Shift-'
        	tkAccel <- sub("shift[+]", "Shift-", tkAccel)
        	## 'alt+' becomes 'Alt-'
        	tkAccel <- sub("alt[+]", "Alt-", tkAccel)
			## Get parent window name
			pWin <- sub("^[$]Tk[.]([a-zA-Z0-9 _.-]+)/.*$", "\\1", menu)
			## Create the binding
			cmd <- paste('tkbind(WinGet("', pWin, '"), "', tkAccel,
				'", function() tkMenuInvoke("', menu, '", "', item, '"))',
				sep = "")
			eval(parse(text = cmd))
		} else Aopt <- ""
		## Rework options
		if (options == "") opts <- "" else opts <- paste(",", options)
		cmd <- paste('tkadd(M, "command", label = "', lbl,
			'", command = function() tkMenuItemCall( { ', action,
			' } ) , compound = "left"', Iopt, Aopt, Uopt, opts, ')', sep = "")
		eval(parse(text = cmd))
	}
	## Register this menu entry in .guiMenus
	entry <- data.frame(name = I(item), action = I(action), image = I(image),
		accel = I(accel), options = I(options))
	items <- .guiMenus[[menu]]$Items
	if (is.null(items)) {
		.guiMenus[[menu]]$Items <- entry
	} else {
        .guiMenus[[menu]]$Items <- rbind(items, entry)
	}
	## Update the SciViews:TempEnv version of .guiMenus
	assignTemp(".guiMenus", .guiMenus)
	invisible(item)
}

tkMenuDelItem <- function (menu, item)
{
	## Look if the menu exists (must be in .guiMenus in SciViews:TempEnv)
	.guiMenus <- getTemp(".guiMenus")
	M <- .guiMenus[[menu]]
	if (is.null(M))
		return(invisible(FALSE))
	## Look if the item exists
	Items <- M$Items
	Pos <- which(Items$name == item) - 1  # Because Tk menu indices start at 0
	if (length(Pos) == 0)
		return(invisible())
	## Check that this item is not a submenu (must use tkMenuDel() instead)
	if (Items$action[Pos + 1] == "[menu]")
		stop("item ", item, " is a submenu. Use tkMenuDel() instead!")
	## Delete that entry
	tkdelete(M, Pos)
	## Eliminate that entry from .guiMenus
	Items <- Items[Items$name != item, ]
	.guiMenus[[menu]]$Items <- Items
	## Update the SciViews:TempEnv version of .guiMenus
	assignTemp(".guiMenus", .guiMenus)
	invisible(TRUE)
}

tkMenuDel <- function (menu)
{
	## Delete a whole menu and all submenus
	## Look if the menu exists (must be in .gui.Menus in SciViews:TempEnv)
	.guiMenus <- getTemp(".guiMenus")
	M <- .guiMenus[[menu]]
	if (is.null(M))
		return(invisible(FALSE))
	## Look at all submenus to delete
	Menus <- names(.guiMenus)
	Mmatch <- (substr(Menus, 1, nchar(menu)) == menu)
	dMenus <- sort(Menus[Mmatch], decreasing = TRUE)  # Sort menus bottom to top
	## Delete each menu in turn
	for (i in 1:length(dMenus)) {
		Pname <- dirname(dMenus[i])
		P <- .guiMenus[[Pname]]
		N <- basename(dMenus[i])
		Items <- P$Items
		Pos <- which(P$Items$name == N)
		## If this is not a toplevel menu, Pos is Pos - 1
### TODO: consider also tearoff menus that way!
		if (regexpr("/", Pname) > 0) Pos <- Pos - 1
		if (length(Pos) > 0) {
			tkdelete(P, Pos)
			.guiMenus[[Pname]]$Items <- Items[Items$name != N, ]
		}
		.guiMenus[[dMenus[i]]] <- NULL  # Eliminate the entry
	}
	## Update the SciViews:TempEnv version of .guiMenus
	assignTemp(".guiMenus", .guiMenus)
	invisible(TRUE)
}

tkMenuChangeItem <- function (menu, item, action = "", options = "")
{
	## The Tk version of MenuChangeItem()
	## Look if the menu exists (must be in .guiMenus in SciViews:TempEnv)
	.guiMenus <- getTemp(".guiMenus")
	M <- .guiMenus[[menu]]
	if (is.null(M))
		return(invisible(FALSE))
	## Look if the item exists
	Items <- M$Items
	Pos <- which(Items$name == item) - 1  # Because Tk menu indices start at 0
	if (length(Pos) == 0)
		return(invisible())
	## Check that this item is not a submenu or a separator
	Act <- Items$action[Pos + 1]
	if (Act == "[separator]")
		stop("item ", item, " is a separator; its state cannot be changed!")
    if (Act == "[menu]")
		stop("item ", item, " is a menu; its state cannot be changed!")
	if (options != "") {
		## Change the configuration of that entry
		eval(parse(text = paste("tkentryconfigure(M, Pos, ", options, ")",
			sep = "")))
	}
	## Do we need to change the action?
	if (action != "" && action != Act) {
		## Change the action
        cmd <- paste('tkentryconfigure(M, Pos, command = function() ',
			action, ')', sep = "")
		eval(parse(text = cmd))
		## Update .guiMenus
		Items$action[Pos + 1] <- action
		.guiMenus[[menu]]$Items <- Items
		## Update the SciViews:TempEnv version of .guiMenus
		assignTemp(".guiMenus", .guiMenus)
		return(invisible(TRUE))
	}
}

tkMenuStateItem <- function (menu, item, active = TRUE)
{
	## The Tk version of MenuStateItem()
	## Look if the menu exists (must be in .guiMenus in SciViews:TempEnv)
	.guiMenus <- getTemp(".guiMenus")
	M <- .guiMenus[[menu]]
	if (is.null(M))
		return(invisible(FALSE))  # If menu does not exists!
	## Look if the item exists
	Items <- M$Items
	Pos <- which(Items$name == item) - 1  # Because Tk menu indices start at 0
	if (length(Pos) == 0)
		return(invisible())
	## Check that this item is not a separator
	if (Items$action[Pos + 1] == "[separator]")
		stop("item ", item, " is a separator; its state cannot be changed!")
	## Set state for that entry
	State <- if (active) "normal" else "disabled"
	tkentryconfigure(M, Pos, state = State)
	invisible(active)
}

tkMenuInvoke <- function (menu, item)
{
	## Given a menu and an item in this menu, trigger the item action
	## Look if the menu exists (must be in .guiMenus in SciViews:TempEnv)
	.guiMenus <- getTemp(".guiMenus")
	M <- .guiMenus[[menu]]
	if (is.null(M))
		return(invisible(FALSE))
	## Look if the item exists
	Items <- M$Items
	Pos <- which(Items$name == item) - 1  # Because Tk menu indices start at 0
	if (length(Pos) == 0)
		return(invisible())
	## Invoke this menu entry
	tcl(M, "invoke", Pos)
	invisible(TRUE)
}
