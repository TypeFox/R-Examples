sendSocketClients <- function (text, sockets = "all", serverport = 8888)
{
    ## Note that 'real' clients should manage to print this BEFORE the current
    ## command line, something that 'SimpleClient.Tcl' cannot do!

    ## Make sure that the text ends with a carriage return
    ## (same behavior as in Mac R.app but different from RGui!)
    if (regexpr("\n^", text) < 0) text <- paste(text, "\n", sep = "")

    ## Send the given text to one or more clients through a socket
    if (sockets == "all")
		sockets <- getSocketClientsNames(port = serverport)
    if (!is.null(sockets) && length(sockets) > 0)
		for (i in 1:length(sockets))
			tcl("puts", sockets[i], text)
}
