winMenuChangeItem <- function (menu, item, action, options = "")
{
### TODO: this is buggy under R 2.2.0! Recheck for latest R version
	if (action == "")
		return(invisible())
	if (!isRgui())
		stop("This function can only be used with Rgui under Windows!")
	## First check if the entry exists
	if (!menu %in% winMenuNames())
		stop("menu '", menu, "' does not exist!")
	if (!item %in%  names(winMenuItems(menu)))
		stop("item '", item, "' does not exist!")
	if (action == "enable" || action == "disable")
		stop("Use 'winMenuStateItem()' instead ",
			"to enable/disable winMenu entries!")
	winMenuAddItem(menu, item, action)
	if (options != "") {
		res <- switch(options,
			"enable" = winMenuAddItem(menu, item, "enable"),
			"disable" = winMenuAddItem(menu, item, "disable"),
			warning("'options = ", options, "' is not supported")
		)
	}
	invisible(TRUE)
}

winMenuStateItem <- function (menu, item, active = TRUE)
{
    if (!isRgui())
		stop("This function can only be used with Rgui under Windows!")
	## First check if the entry exists
	if (!menu %in% winMenuNames())
		stop("menu '", menu, "' does not exist!")
	if (!item %in%  names(winMenuItems(menu)))
		stop("item '", item, "' does not exist!")
	## Enable/disable the entry
	cmd <- if (active) "enable" else "disable"
	winMenuAddItem(menu, item, cmd)
	invisible(TRUE)
}

winMenuInvoke <- function (menu, item)
{
	if (!menu %in% winMenuNames())
		stop("menu '", menu, "' does not exist!")
	action <- winMenuItems(menu)[item]
	if (is.null(action) || is.na(action) || length(action) == 0 || action == "")
		return(invisible(FALSE))
	## Print the action on the console
	cat(getOption("prompt"), action, "\n", sep = "")
	## Add the action to the command history
	timestamp(action, prefix = "", suffix = "", quiet = TRUE)
	## and execute it in .GlobalEnv
	eval(parse(text = action), envir = .GlobalEnv)
	invisible(TRUE)
}
