tkToolItems <- function (toolbar)
{
	T <- getTemp(".guiTools")[[toolbar]]
	if (is.null(T))
		stop("unable to retrieve items for ", toolbar, "\n(toolbar ", toolbar,
			" does not exist)")
	Titems <- T$Items
	if (is.null(Titems)) {
		res <- character(0)
	} else {
		res <- as.character(Titems$action)
	}
	if (length(res) > 0) names(res) <- Titems$name
	res
}

tkToolAdd <- function (toolbar, side = "top")
{
	## Get the toolbar name
	Tname <- basename(toolbar)
	## Get the name of the parent
	Pname <- dirname(toolbar)
	## Look if the parent exists (must be in .guiTools in SciViews:TempEnv)
	.guiTools <- getTemp(".guiTools")
	## Do not create the toolbar if it already exists
	if (toolbar %in% names(.guiTools))
		stop("toolbar ", toolbar, " already exists!")
	Parent <- .guiTools[[Pname]]
	if (is.null(Parent)) {
        ## If base toolbar is a "root", try to create the corresponding
		## Tk toolbar area
        if (regexpr("/", Pname) < 0) {
			## Get the name of the Tk window that hosts the menu
			TkWinName <- sub("^[$]Tk[.]", "", Pname)
			## Look if such a Tk window is recorded in .guiWins
			TkWin <- winGet(TkWinName)
			if (is.null(TkWin))
				stop("toolbar does not exist, ",
					"and the parent Tk window cannot be found!")
			## Create the toolbar
            Parent <- ttkframe(TkWin)
			tkpack(Parent, side = side, anchor = "w")
			## Add an entry for this top menu in .guiTools
			.guiTools[[Pname]] <- Parent
		} else stop("unable to add toolbar\n(base toolbar does not exists)")
	}
	## Check that the name is not already used (e.g., for a toolbar entry)
    items <- Parent$Items
    if (!is.null(items) && Tname %in% items$name)
		stop("the toolbar name is already used in this toolbar!")
	## Now that the toolbar area exists, we can create this toolbar
	item <- Tname
	Child <- ttkframe(Parent)
    tkpack(Child, side = side, anchor = "w")
 	## Add an entry for this toolbar menu in .guiTools
	.guiTools[[toolbar]] <- Child
    ## ... and register it in the items of parent toolbars in .guiTools
	entry <- data.frame(name = I(Tname), action = I("[toolbar]"), image = I(""),
		side = I(side))
	if (is.null(items)) {
		.guiTools[[Pname]]$Items <- entry
	} else {
        .guiTools[[Pname]]$Items <- rbind(items, entry)
	}
	## Update the SciViews:TempEnv version of .guiTools
	assignTemp(".guiTools", .guiTools)
	invisible(toolbar)
}

tkToolAddItem <- function (toolbar, item, action, image = "", options = "")
{
	## Look if the toolbar exists (must be in .guiTools in Sciviews:TempEnv)
	.guiTools <- getTemp(".guiTools")
	Tl <- .guiTools[[toolbar]]
	if (is.null(Tl))
		stop("toolbar does not exist!")
	
	## Is this item a separator?
	isSeparator <- regexpr("^-+$", item) > 0
	
	Items <- Tl$Items
	## Look if the toolbar already has this tool
	if (is.null(Items)) n <- 0 else {
		if (!isSeparator && item %in% Items$name)
			stop("the tool ", item, " already exists!")
		n <- max(Items$pos) + 1
	}
		
	## Add the entry at the end of the Tk toolbar (### TODO: allow other positions)
	## First look if it is a tool button or a separator
	if ( isSeparator ) {  # This must be a separator
### TODO: choose correct orientation
		but <- ttkseparator(Tl, orient = "vertical")
		tkgrid(but, row = 0, column = n, sticky = "nsew")
		action <- "[separator]"
		options <- 'orient = "vertical"'
	} else {  # This is a tool button	
		## Rework options
### TODO: I need to deal with options!
		##if (options == "") opts <- "" else opts <- paste(",", options)
		## Do we have to add an image to the tool button?
		if (image != "") {
			Img <- imgGet(image)
			but <- ttkbutton(Tl, text = item, image = as.character(Img),
				compound = "image", style = "Toolbutton",
				command = .actionWrapper(action))
				# command = eval(parse(text = paste("function()",  action))))
		} else {
			but <- ttkbutton(Tl, text = item,
				compound = "left", style = "Toolbutton",
				command = .actionWrapper(action)) 
				# command = eval(parse(text = paste("function()",  action))))
		}
		tkgrid(but, row = 0, column = n, sticky = "nsew")
### TODO: This needs tcltk2 => how to get rid of this dependency?
		#tk2tip(but, lbl)
	}
	## Add an entry for this toolbutton in .guiTools
	itempath <- paste(toolbar, item, sep = "/")
	.guiTools[[itempath]] <- but
	## Register this tool in .guiTools
	entry <- data.frame(name = I(item), action = I(action), image = I(image),
		options = I(options), pos = n)
	items <- .guiTools[[toolbar]]$Items
	if (is.null(items)) {
		.guiTools[[toolbar]]$Items <- entry
	} else {
        .guiTools[[toolbar]]$Items <- rbind(items, entry)
	}
	## Update the SciViews:TempEnv version of .guiTools
	assignTemp(".guiTools", .guiTools)
	invisible(itempath)
}

.actionWrapper <- function (action)
{
	function () {
		timestamp(action, prefix = "", suffix = "") # Was: .Internal(addhistory(action))
		cat(">> ", action, "\n") 
		eval(parse(text = action))
	}
}


tkToolDelItem <- function (toolbar, item)
{
	## Look if the toolbar exists (must be in .guiTools in SciViews:TempEnv)
	.guiTools <- getTemp(".guiTools")
	Tl <- .guiTools[[toolbar]]
	if (is.null(Tl)) return(invisible(FALSE))
	## Look if the item exists
	Items <- Tl$Items
	Pos <- which(Items$name == item)
	if (length(Pos) == 0) return(invisible(FALSE))
	## Check that this item is not a submenu (must use tkMenuDel() instead)
	if (Items$action[Pos] == "[menu]")
		stop("item ", item, " is a submenu. Use tkMenuDel() instead!")
	## Delete that widget
	itempath <- paste(toolbar, item, sep = "/")
	tkdestroy(.guiTools[[itempath]])
	.guiTools[[itempath]] <- NULL
	## Eliminate that entry from .guiTools
	Items <- Items[-Pos, ]
	.guiTools[[toolbar]]$Items <- Items
	## Update the SciViews:TempEnv version of .guiTools
	assignTemp(".guiTools", .guiTools)
	return(invisible(TRUE))
}

tkToolDel <- function (toolbar)
{
	## Delete a whole toolbar
	## Look if the toolbar exists (must be in .guiTools in Sciviews:TempEnv)
	.guiTools <- getTemp(".guiTools")
	Tl <- .guiTools[[toolbar]]
	if (is.null(Tl)) return(invisible())
	## Look at all toolbuttons to delete
	Toolbars <- names(.guiTools)
	Tmatch <- (substr(Toolbars, 1, nchar(toolbar)) == toolbar)
	dTools <- sort(Toolbars[Tmatch], decreasing = TRUE) # Sort toolbars decr.
	## Delete each tool
	for (i in 1:length(dTools)) {
		tkdestroy(.guiTools[[dTools[i]]])
		.guiTools[[dTools[i]]] <- NULL     # Eliminate the entry
	}
	## Eliminate the entry for the toolbar in the items list
	Items <- .guiTools[[dirname(toolbar)]]$Items
	Pos <- which(Items$name == basename(toolbar))
	if (length(Pos) > 0) Items <- Items[-Pos, ]
	.guiTools[[dirname(toolbar)]]$Items <- Items
	## Update the SciViews:TempEnv version of .guiTools
	assignTemp(".guiTools", .guiTools)
	return(invisible(TRUE))
}

tkToolChangeItem <- function (toolbar, item, action = "", options = "")
{
	## The Tk version of ToolChangeItem()
	## Look if the toolbar exists (must be in .guiTools in SciViews:TempEnv)
	.guiTools <- getTemp(".guiTools")
	Tl <- .guiTools[[toolbar]]
	if (is.null(Tl)) return(invisible(FALSE))
	## Look if the item exists
	Items <- Tl$Items
	itempath <- paste(toolbar, item, sep = "/")
	Pos <- which(Items$name == item)
	if (length(Pos) == 0) return(invisible())
	## Check that this item is not a submenu or a separator
	Act <- Items$action[Pos]
	if (Act == "[separator]")
		stop("item ", item, " is a separator; its state cannot be changed!")
    if (Act == "[menu]")
		stop("item ", item, " is a menu; its state cannot be changed!")
	if (options != "") {
		## Change the configuration of that entry
		eval(parse(text = paste('tkconfigure("', as.character(.guiTools[[itempath]]),
			'", Pos, ', options, ')', sep = "")))
		Items$options[Pos] <- options
	}
	## Do we need to change the action?
	if (action != "" && action != Act) {
		# Change the action
        tkconfigure(.guiTools[[itempath]], command = eval(parse(text = paste("function()",
			action))))
		## Update .guiTools
		Items$action[Pos] <- action
	}
	.guiTools[[toolbar]]$Items <- Items
	## Update the SciViews:TempEnv version of .guiTools
	assignTemp(".guiTools", .guiTools)
	return(invisible(TRUE))
}

tkToolStateItem <- function (toolbar, item, active = TRUE)
{
	## The Tk version of ToolStateItem()
	## Look if the toolbar exists (must be in .guiTools in SciViews:TempEnv)
	.guiTools <- getTemp(".guiTools")
	Tl <- .guiTools[[toolbar]]
	if (is.null(Tl)) return(invisible())
	## Look if the item exists
	Items <- Tl$Items
	itempath <- paste(toolbar, item, sep = "/")
	Pos <- which(Items$name == item)
	if (length(Pos) == 0) return(invisible(FALSE))
	## Check that this item is not a separator
	if (Items$action[Pos] == "[separator]")
		stop("item ", item, " is a separator; its state cannot be changed!")
	## Set state for that entry
	State <- if (active) "normal" else "disabled"
	tkconfigure(.guiTools[[itempath]], state = State)
	return(invisible(active))
}

tkToolInvoke <- function (toolbar, item)
{
	## Given a toolbar and an item in this toolbar, trigger the item action
	## Look if the toolbar exists (must be in .guiTools in SciViews:TempEnv)
	.guiTools <- getTemp(".guiTools")
	Tl <- .guiTools[[toolbar]]
	if (is.null(Tl)) return(invisible(FALSE))
	## Look if the item exists
	Items <- Tl$Items
	itempath <- paste(toolbar, item, sep = "/")
	Pos <- which(Items$name == item)
	if (length(Pos) == 0) return(invisible())
	## Invoke this tool entry
	tkinvoke(.guiTools[[itempath]])
	return(invisible(TRUE))
}
