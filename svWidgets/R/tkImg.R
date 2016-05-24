tkImgAdd <- function (file, type = "gif", update = FALSE)
{
	## Add a Tk image to the list (just GIF for the moment,
	## but the Tcl/Tk Img package allows for more!)
	if (type != "gif")
		stop("Only 'gif' images currently supported!")
	if (!file.exists(file))
		stop("File '", file, "' not found!")
	## Load the image and assign it to an item of .guiImgs in SciViews:TempEnv
	.guiImgs <- getTemp(".guiImgs")
	if (is.null(.guiImgs)) {
		.guiImgs <- list()
		class(.guiImgs) <- c("guiImg", "gui", class(.guiImgs))
	}
	## Calculate image name as being the basename of the file without extension
	Iname <- sub("[.][^.]+$", "", basename(file))
	## In order to indicate it is a Tk resource, prepend '$Tk.'
	Iname <- paste("$Tk.", Iname, sep = "")
	## If that name already exists, do nothing, except if we ask to update it
	if (Iname %in% names(.guiImgs)) {
		if (isTRUE(update)) {
			## Delete the image to recreate it with new resource
			tcl("image", "delete", Iname)
		} else return(invisible(Iname))	# Do nothing
	}
	Image <- tclVar()
	tcl("image", "create", "photo", Image, file = file)

	.guiImgs[[Iname]] <- Image
	## Reassign .guiImgs to SciViews:TempEnv
	assignTemp(".guiImgs", .guiImgs)
	invisible(Iname)
}

tkImgDel <- function (image)
{
	## Delete a tk image ressource from the list
	.guiImgs <- getTemp(".guiImgs")
	## Is the image there?
	if (!image %in% names(.guiImgs))
		return(invisible(FALSE))
	## Delete the image
	Image <- .guiImgs[[image]]
	tcl("image", "delete", Image)
	## Eliminate it from the list in .guiImgs
	.guiImgs[[image]] <- NULL
	## Reassign .guiImgs to SciViews:TempEnv
	assignTemp(".guiImgs", .guiImgs)
	## Indicate that the image is actually deleted
	invisible(TRUE)
}

tkImgRead <- function (dir, type = "gif")
{
	## Read all gif images from a directory into tkImage resources
	## Check that the dir exists
	if (!file.exists(dir) || !file.info(dir)$isdir)
		stop("'", dir, "' does not exist, or is not a directory!")
	## List all file of 'type' in that directory
	if (type != "gif")
		stop("only type = 'gif' is currently supported")
	pattern <- "[.][gG][iI][fF]$"
	files <- list.files(dir, pattern = pattern, full.names = TRUE)
	if (length(files) == 0)
		return(invisible())
	for (i in 1:length(files))
		tkImgAdd(files[i], type = type)
	invisible(files)
}
