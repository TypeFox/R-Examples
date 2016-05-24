print.guiTool <- function (x, ...)
{
	cat("A SciViews GUI tool object:", "\n")
	print(unclass(x))
	invisible(x)
}

toolAdd <- function (toolbar, side = "top")
{
	invisible(switch(toolType(toolbar),
		tkTool = tkToolAdd(toolbar = toolbar, side = side)))
}

toolAddItem <- function (toolbar, item, action, image = "", options = "")
{
	invisible(switch(toolType(toolbar),
		tkTool = tkToolAddItem(toolbar = toolbar, item = item, action = action,
			image = image, options = options)))
}

toolDel <- function (toolbar)
{
    invisible(switch(toolType(toolbar),
		tkTool = tkToolDel(toolbar = toolbar)))
}

toolDelItem <- function (toolbar, item)
{
    invisible(switch(toolType(toolbar),
		tkTool = tkToolDelItem(toolbar = toolbar, item = item)))
}

toolNames <- function ()
{
	res <- character(0)
	## Retrieve toolbar names from tk toolbars
	res <- c(res, names(getTemp(".guiTools")))
	## Eliminate toplevel entries
    if (length(res) > 0) res <- res[regexpr("/", res) > 0]
	res
}

toolItems <- function (toolbar)
{
    switch(toolType(toolbar),
		tkTool = tkToolItems(toolbar = toolbar))
}

toolType <- function (toolbar, warn = TRUE)
{
	## Given a toolbar, return its type ("tkTool", NA)
	if (regexpr("^[$]Tk[.].+/", toolbar) > 0) {
			res <- "tkTool"
	} else {
        if (warn) warning("Unrecognized toolbar type for ", toolbar)
		res <- NA
	}
	res
}

toolChangeItem <- function (toolbar, item, action = "", options = "")
{
	## Change action or options for toolbar entries
	invisible(switch(toolType(toolbar),
		tkTool = tkToolChangeItem(toolbar, item, action, options)))
}

toolStateItem <- function (toolbar, item, active = TRUE)
{
	## Activate/inactivate toolbar entries
	invisible(switch(toolType(toolbar),
		tkTool = tkToolStateItem(toolbar, item, active)))
}

toolInvoke <- function (toolbar, item)
{
	## Invoke a toolbutton
	invisible(switch(toolType(toolbar),
		tkTool = tkToolInvoke(toolbar, item)))
}

toolRead <- function (file = "Tools.txt")
{
	## Read toolbars from a file
    T <- scan(file, character(0), sep = "\n", comment.char = "#", quiet = TRUE)
	## Split the lines into item, command, options
    T <- strsplit(T, "~~")
	## Strip leading and trailing spaces/tabs (used to align items in the file)
    T <- lapply(T, function(x) sub("[ \t]+$", "", sub("^[ \t]+", "", x)))
	## Move line after line and replace '|' by values
    N <- length(T)
	## Must have at least two entries
	if (N < 2) return(invisible())
	## First entry must be a toplevel, thus, it must start with '$'
	if (regexpr("^[$]", T[[1]][1]) < 0)
		stop("first entry is not a toolbar!")
    toolLevels <- T[[1]][1]
	## Initialize a data frame to contain decrypted info
	dat <- rep("", N)
	L <- data.frame(tool = I(dat), item = I(dat), image = I(dat),
		action = I(dat), options = I(dat))
	Litem <- data.frame(tool = I(toolLevels), item = I(""), image = I(""),
		action = I("[toolbar]"), options = I(""))
	L[1, ] <- Litem
	for (i in 2:N) {
		entry <- T[[i]][1]
		## Split on '|'
		split <- strsplit(entry, "[|]")[[1]]
		## Combine toolLevels & split
		last <- length(split)
		toolLevels[last] <- split[last]
		toolLevels <- toolLevels[1:last]
		## Recombine toolLevels for getting recalculated entry
		entry <- paste(toolLevels, collapse = "/")
		## Is this just a tool button, or a menu tool button/menu item?
		lastentry <- basename(entry)
		if (regexpr("^[$]", lastentry) > 0) { # This is a tool button
			## Remove '$' before tool button entries
			tool <- gsub("/[$]", "/", entry)
			item <- ""      	# No item
			image <- ""     	# No image
			action <- if (last == 1) "[toolbar]" else "[tool]"
			options <- ""       # No options (currently)
		} else {  # This is an menu entry in a tool button menu
			## tool = entry minus last item (remove '$' before tool entries)
			tool <- gsub("/[$]", "/", dirname(entry))
			## Decrypt lastentry to get image & item ([image]item)
   			item <- sub("^[[][a-zA-Z0-9 ._-]+[]]", "", lastentry)
			if (item == lastentry) image <- "" else {
				image <- sub("^[[]([a-zA-Z0-9 ._-]+)[]].+$", "\\1", lastentry)
				## Since these are Tk images resources, I have to prefix '$Tk.'
				image <- paste("$Tk.", image, sep = "")
			}
   			action <- T[[i]][2]
			if (is.na(action)) action <- ""
			options <- T[[i]][3]
            if (is.na(options)) options <- ""
		}
        Litem <- data.frame(tool = I(tool), item = I(item), image = I(image),
			action = I(action), options = I(options))
		## Add it to the data.frame
		L[i, ] <- Litem
	}
	
	## The [toolbar] entries are not needed
	L <- subset( L, action != "[toolbar]" )

	## Execute this line-by-line to create the various tools 
	N <- nrow(L)
	for (i in 1:N) {
		action <- L$action[i]
		if (action == "[tool]") { 	# Create a toolbar
            toolAdd(toolbar = L$tool[i])
#		} else if (action == "[tool]") {  # Create a tool button
#			### TODO: determine which type of tool button it is!
#			ToolAddItem(toolbar = L$menu[i])
#		} else {    # Create a menu entry for a menu tool button
#			MenuAddItem(menu = L$tool[i], item = L$item[i], action = L$action[i],
#				image = L$image[i], options = L$options[i])
#		}
		} else {  # Create a tool in the toolbar
            toolAddItem(toolbar = L$tool[i], item = L$item[i],
				action = L$action[i], image = L$image[i],
				options = L$options[i])
		}
	}

	invisible(L)
}

toolReadPackage <- function (package, subdir = "gui", file = "Tools.txt")
{
	## Create toolbars using a toolbar definition file located in a R package
	dir <- system.file(subdir, package = package)
	## Check that the dir exists
	if (!file.exists(dir) || !file.info(dir)$isdir)
		stop("'", dir, "' does not exist, or is not a directory!")
	## Check that the file exists
	File <- file.path(dir, file)
	if (!file.exists(File))
		stop("'", file, "' not found in '", dir, "'!")
	## Read the toolbar file
	invisible(toolRead(File))
}
