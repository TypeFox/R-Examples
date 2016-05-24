koInstall <- function ()
{
## PhG: this is a mechanisms we don't implement at the end!
#	assignTemp(".guiCmd", function(command, ...) {
#		## TODO: define these commands
#		command <- switch(command,
#			load = "",
#			source = "",
#			save = "",
#			import = "",
#			export = "",
#			report = "",
#			setwd = "",
#			"")
#	})
#	# FIXME: these sv.* functions do not exist in SciViews-K!
#	assignTemp(".guiObjBrowse", function(id, data) {
#		koCmd('sv.objBrowse("<<<id>>>", "<<<dat>>>");',
#		list(id = id, dat = data))
#	})
#	assignTemp(".guiObjInfo", function(id, data) {
#		koCmd('sv.objInfo("<<<id>>>", "<<<dat>>>");',
#		list(id = id, dat = data))
#	})
#	assignTemp(".guiObjMenu", function(id, data) {
#		koCmd('sv.objMenu("<<<id>>>", "<<<dat>>>");',
#		list(id = id, dat = data))
#	})
#
#	## Functions specific to Komodo as a GUI client
#	assignTemp(".koCmd", function(command, ...) {
#		## This mechanism avoids dependence on svGUI for packages that provide
#		## functionalities that work with or without Komodo (like svUnit)
#		## Instead of calling koCmd() directly, we look if .koCmd is defined
#		## in tempenv and we run it.
#		## This allows also redefining koCmd() without changing code in the
#		## packages that depend on .koCmd()
#		koCmd(command, ...)
#	})

	## Register a TaskCallback to generate automatically informations for an
	## object browser.
	h <- getTemp(".svTaskCallbackManager", default = NULL, mode = "list")
	if (!is.null(h))
		h$add(koAutoRefresh, name = "koAutoRefresh")
}
