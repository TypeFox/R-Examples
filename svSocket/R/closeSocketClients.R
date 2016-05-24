closeSocketClients <- function (sockets = "all", serverport = 8888)
{
    ## Nicely close socket client(s) by sending "\f"
    ## To be interpreted by a compatible client that manages to close connection
    if (sockets == "all")
		sockets <- getSocketClientsNames(port = serverport)
    if (!is.null(sockets) && length(sockets) > 0)
		for (i in 1:length(sockets))
			tcl("puts", sockets[i], "\n\f")
}
