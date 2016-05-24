print.guiImg <- function (x, ...)
{
	cat("A SciViews GUI image object:", "\n")
	print(unclass(x))
	invisible(x)
}

imgAdd <- function (file, type = "gif", img.type = "tkImage",
update = FALSE, ...)
{
	invisible(switch(img.type,
		tkImage = tkImgAdd(file = file, type = type, update = update),
		stop("Unrecognized image type '", img.type, "'")))
}

imgDel <- function (image)
{
    invisible(switch(imgType(image),
		tkImage = tkImgDel(image)))
}

imgGet <- function (image)
{
	## Get the image
    getTemp(".guiImgs")[[image]]
}

imgType <- function (image, warn = TRUE)
{
	## Get the type of image
	if (regexpr("^[$]Tk[.]", image) > 0) {
		res <- "tkImage"
	} else {
		if (warn) warning("Unrecognized image type for ", image)
		res <- NA
	}
	res
}

imgNames <- function ()
{
	## List all available images
    res <- names(getTemp(".guiImgs"))
	if (is.null(res)) res <- character(0)
	res
}

imgRead <- function (dir, type = "gif", img.type = "tkImage")
{
	## Depending on 'img.type', we call a different function
	invisible(switch(img.type,
		tkImage = tkImgRead(dir = dir, type = type),
		stop("Unrecognized image type '", img.type, "'")))
}

imgReadPackage <- function (package, subdir = "gui", type = "gif",
img.type = "tkImage")
{
	## Create image resources by reading a series of image files from the 'gui'
	## subdirectory of a package
	dir <- system.file(subdir, package = package)
	invisible(imgRead(dir = dir, type = type, img.type = img.type))
}
