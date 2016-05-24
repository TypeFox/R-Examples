koUninstall <- function ()
{
## PhG: this is a mechanisms we don't implement at the end...
#	## Eliminate .guiCmd
#	rmTemp(".guiCmd")
#	rmTemp(".guiObjBrowse")
#	rmTemp(".guiObjInfo")
#	rmTemp(".guiObjMenu")
#
#	rmTemp(".koCmd")

	## Unregister our own TaskCallback
	h <- getTemp(".svTaskCallbackManager", default = NULL, mode = "list")
	if (!is.null(h))
		try(h$remove("koAutoRefresh"), silent = TRUE)
}
