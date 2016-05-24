.onLoad <- function (lib, pkg)
{
	## TODO: check configuration and install everything that we need to use the
	## SciViews extensions, including the HTTP or socket server
	#serve <- getOption("ko.serve")
	#if (!is.null(serve)) {
	#	startSocketServer(port = as.integer(serve)[1])
	#	guiInstall()
	#}
}

.onUnload <- function (libpath)
{
	#serve <- getOption("ko.serve")
	#if (!is.null(serve) && "package:svSocket" %in% search())
	#	stopSocketServer(port = as.integer(serve)[1])
	#guiUninstall()
}

.packageName <- "SciViews"
