getSocketServers <- function ()
{
    ## Get the list of currently running socket servers
    return(TempEnv()$SocketServers)
}

getSocketClients <- function (port = 8888)
{
    if (!is.numeric(port[1]) || port[1] < 1)
        stop("'port' must be a positive integer!")
    portnum <- round(port[1])
    port <- as.character(portnum)

    ## Does a server exist on this port?
    servers <- getSocketServers()
    if (!(port %in% servers))
		return(NULL)  # If no R socket server running on this port

    ## Get the list of clients currently connected to this server
    clients <- as.character(.Tcl(paste("array names Rserver", port, sep = "_")))
    ## Eliminate "main", which is the connection socket
    clients <- clients[clients != "main"]

    ## Are there client connected?
    if (length(clients) == 0) return(character(0))

    ## For each client, retrieve its address and port
    addresses <- NULL
    arrayname <- paste("Rserver", port, sep = "_")
    for (i in 1:length(clients)) {
		client <- as.character(.Tcl(paste("array get", arrayname, clients[i])))
		addresses[i] <- sub(", ", ":", client[2])
    }
    names(addresses) <- clients
    return(addresses)
}

getSocketClientsNames <- function (port = 8888)
    names(getSocketClients(port = port))

getSocketServerName <- function (port = 8888)
{
    if (!is.numeric(port[1]) || port[1] < 1)
        stop("'port' must be a positive integer!")
    portnum <- round(port[1])
    port <- as.character(portnum)

    ## Return the name of a given R socket server
    servers <- getSocketServers()
    if (!(port %in% servers))
		return(NULL)  # If no R socket server running on this port

    ServerNames <- names(servers)
    return(ServerNames[servers == port])
}
