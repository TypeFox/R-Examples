.onLoad <- function (lib, pkg)
{
	## Create our SciViews task callback manager
	#PhG: now moved to svKomodo!
	#assignTemp(".svTaskCallbackManager", svTaskCallbackManager())
}

.onUnload <- function (libpath)
{
	#PhG: now moved to svKomodo!
	#removeTaskCallback("SV-taskCallbackManager")
	#rmTemp(".svTaskCallbackManager")
	
	#PhG: From .Last.lib(), now in .onUnload()
	## Make sure that all clients are disconnected
    ## and all servers are closed
    Servers <- getSocketServers()
    if (is.null(Servers) || length(Servers) < 1) return()
    cat(ngettext(length(Servers), "Stopping socket server\n",
        "Stopping socket servers\n"))
    stopSocketServer("all")
    ## TODO: make sure to delete all client environments
	## (or do it in stopSocketServer()?)
}
