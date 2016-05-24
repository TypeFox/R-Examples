stopSocketServer <- function (port = 8888)
{
    ## Stop one or more running socket server(s)
    if (port == "all") {
		port <- getSocketServers()
		servers <- port
    } else servers <- getSocketServers()
    if (!is.numeric(port) || any(port < 1))
        stop("'port' must be positive integers!")
    port <- round(port)
    anyclosed <- FALSE
    for (i in 1:length(port)) {
		Port <- port[i]
		if (Port %in% servers) {  # This port is open
			anyclosed <- TRUE
			## First ask to all clients to nicely disconnect (note: if they don't
			## the server simply does not process them any more!)
			closeSocketClients(serverport = Port)

			## Assign it back, with the corresponding port stripped out
			## But if I was the last one, delete the SocketServers variable
			servers <- servers[servers != Port]
			if (length(servers) == 0) {
				if (exists("SocketServers", envir = TempEnv(),
					inherits = FALSE)) rm("SocketServers", envir = TempEnv())
			} else {
				assign("SocketServers", servers[servers != Port],
					envir = TempEnv())
			}

			## Eliminate the processing function from SciViews:TempEnv
			sockProc <- paste("SocketServerProc", Port, sep = "_")
			if (exists(sockProc, envir = TempEnv()))
				rm(list = sockProc, envir = TempEnv())

			## Close the socket in order not to reject future client connections
			.Tcl(paste("close $Rserver_", Port, "(main)", sep = ""))

			## Note: Tcl procs and variables are not eliminated yet
			## because there may be still clients connected!
		}
    }
    return(anyclosed)
}
