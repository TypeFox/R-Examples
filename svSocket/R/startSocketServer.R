startSocketServer <- function (port = 8888, server.name = "Rserver",
procfun = processSocket, secure = FALSE, local = !secure)
{
    ## OK, could be port = 80 to emulate a simple HTML server
    ## This is the main function that starts the server
    ## This function implements a basic R socket server on 'port'
    ## SocketServerProc is the R workhorse function that do the computation
    ## The server is written in Tcl. This way it is not blocking R command-line!
    ## It is designed in a way that R can open simultaneously several ports and
    ## accept connection from multiple clients to each of them.
    ## Commands from each port can be processed differently

	## Secure server requires the tcl-tls package!
	if (isTRUE(secure)) {
		## TODO: On Mac with AquaTclTk installed, I need: addTclPath("/System/Library/Tcl")
		res <- tclRequire("tls")
		if (!inherits(res, "tclObj"))
			stop("You must install the tcl-tls package for using a secure server!")
	}

    if (!is.function(procfun))
		stop("'procfun' must be a function!")
    
	## Note: the data send by the client is in the Tcl $::sockMsg variable
    ## Could a clash happen here if multiple clients send data at the
    ## same time to the R socket server???
    if (!is.numeric(port[1]) || port[1] < 1)
        stop("'port' must be a positive integer!")
    portnum <- round(port[1])
    port <- as.character(portnum)

    if (!is.character(server.name))
		stop("'server.name' must be a string!")
    server.name <- as.character(server.name)[1]

    ## Check if the port is not open yet
    servers <- getSocketServers()
    if (port %in% servers) return(TRUE)  # This port is already open!

    ## We need Tcl to be able to call an R function to process clients' requests
    "tclProcExists" <- function (proc) {
		proc <- as.character(proc[1])
		return(length(as.character(tcl("info", "commands", proc))) == 1)
    }

    if (!tclProcExists("SocketServerProc")) {
		## Create the callback when a client sends data
		"SocketServerFun" <- function () {
			## Note: I don't know how to pass arguments here.
			## So, I use Tcl global variables instead:
			## - the server port from $::sockPort,
			## - the socket client from $::sockClient,
			## - and the message from $::sockMsg
			"tclGetValue_" <- function (name) {
				## Get the value stored in a plain Tcl variable
				if (!is.character(name)) stop("'name' must be a character!")
	
				## Create a temporary dual variable with tclVar()
				Temp <- tclVar(init = "")

				## Copy the content of the var of interest to it
				.Tcl(paste("catch {set ", as.character(Temp), " $", name, "}",
					sep = ""))

				## Get the content of the temporary variable
				Res <- tclvalue(Temp) # Temp is destroyed when function exists
				return(Res)
			}
			
			"TempEnv_" <- function () {
				pos <-  match("SciViews:TempEnv", search())
				if (is.na(pos)) {  # Must create it
					`SciViews:TempEnv` <- list()
					Attach <- function (...) get("attach", mode = "function")(...)
					Attach(`SciViews:TempEnv`, pos = length(search()) - 1)
					rm(`SciViews:TempEnv`)
					pos <- match("SciViews:TempEnv", search())
				}
				return(pos.to.env(pos))
			}

			"getTemp_" <- function (x, default = NULL, mode = "any") {
				if  (exists(x, envir = TempEnv_(), mode = mode,
						inherits = FALSE)) {
					return(get(x, envir = TempEnv_(), mode = mode,
							inherits = FALSE))
				} else {  # Variable not found, return the default value
					return(default)
				}
			}

			"process" <- function () {
				port <- tclGetValue_("::sockPort")
				if (port == "") return(FALSE)  # The server is closed
				client <- tclGetValue_("::sockClient")
				if (client == "") return(FALSE)  # The socket client is unknown!
				msg <- tclGetValue_("::sockMsg")
				if (msg == "") return(FALSE)  # No message!

				## Make sure this message is not processed twice
				.Tcl("set ::sockMsg {}")

				## Do we have to debug socket transactions
				Debug <- isTRUE(getOption("debug.Socket"))
				if (Debug) cat(client, " > ", port, ": ", msg, "\n", sep = "")

				## Function to process the client request: SocketServerProc_<port>
				proc <- getTemp_(paste("SocketServerProc", port, sep = "_"),
					mode = "function")
				if (is.null(proc) || !is.function(proc))
					return(FALSE)  # The server should be closed
				## Call this function
				res <- proc(msg, client, port)
				## Return result to the client
				if (res != "") {
					if (Debug) cat(port, " > ", client, ": ", res, "\n", sep = "")
					chk <- try(tcl("puts", client, res), silent = TRUE)
					if (inherits(chk, "try-error")) {
						warning("Impossible to return results to a disconnected client.")
						return(FALSE)
					}
				}
				return(TRUE)  # The command is processed
			}
			return(process)  # Create the closure function for .Tcl.callback()
		}
		assignTemp("SocketServerProc", SocketServerFun())
		## Create a Tcl proc that calls this function back
		res <- .Tcl.callback(getTemp("SocketServerProc"), TempEnv())
		if (length(grep("R_call ", res) > 0)) {
			## Create a proc with the same name in Tcl
			.Tcl(paste("proc SocketServerProc {} {", res, "}", sep = ""))
		} else stop("Cannot create the SciViews socket server callback function")
    }

    ## Copy procfun into SciViews:TempEnv as SocketServerProc_<port>
    assign(paste("SocketServerProc", port, sep ="_"), procfun, envir = TempEnv())

    ## Create the Tcl function that retrieves data from the socket
    ## (command send by the client), call the processing R function
    ## and returns result to the client
    cmd <- paste(c(paste("proc  sockHandler_", port, " {sock} {", sep = ""),
        paste("global Rserver_", port, sep = ""),
		"if {[eof $sock] == 1 || [catch {gets $sock line}]} {",
		"    # end of file or abnormal connection drop",
		"    fileevent $sock readable {}",
		"    close $sock",
		paste("    #puts \"Close $Rserver_", port, "($sock)\"", sep = ""),
		paste("    unset Rserver_", port, "($sock)", sep = ""),
		"} else {",
		"    # Do we have to redirect the connection?",
		"    if {[string compare \">>>>>>sock\" [string range $line 0 9]] == 0} {",
		"        set redirSock [string range $line 6 12]",
		"        fileevent $sock readable [list sockRedirect $sock $redirSock]",
		paste("        unset Rserver_", port, "($sock)", sep = ""),
		"    } else {",
		"        global sockPort",
		"        global sockClient",
		"        global sockMsg",
		paste("        set ::sockPort", port),
		"        set ::sockClient $sock",
		"        set ::sockMsg $line",
		"        SocketServerProc    ;# process the command in R",
		"}\n}\n}"),
	collapse = "\n")
    ## if {[gets $sock line] < 0} {return} # To handle incomplete lines!
    .Tcl(cmd)

    ## Create the Tcl function that accepts input from a client
    ## (a different one for each server port)
	## Code is slightly different if the server is only local or not
	if (isTRUE(local)) {
		cmd <- paste(c(paste("proc sockAccept_", port, " {sock addr port} {",
			sep = ""),
			paste("global Rserver_", port, sep = ""),
			"# Configure the socket",
			"fconfigure $sock -buffering line -blocking 0",
			"# Accept only local clients",
			"if {$addr != \"127.0.0.1\"} {",
			" #   puts $sock \"Error: Only local clients allowed!\"",
			"    close $sock",
			"    return",
			"}",
			paste("set Rserver_", port, "($sock) [list $addr, $port]", sep = ""),
			paste("fileevent $sock readable [list sockHandler_", port,
				" $sock] }", sep = "")),
		collapse = "\n")
	} else {
		cmd <- paste(c(paste("proc sockAccept_", port, " {sock addr port} {",
			sep = ""),
			paste("global Rserver_", port, sep = ""),
			"# Configure the socket",
			"fconfigure $sock -buffering line -blocking 0",
			paste("set Rserver_", port, "($sock) [list $addr, $port]", sep = ""),
			paste("fileevent $sock readable [list sockHandler_", port,
				" $sock] }", sep = "")),
		collapse = "\n")
	}
	.Tcl(cmd)
	
	## Create a Tcl procedure to redirect output (used in socketClientConnection())
	if (!tclProcExists("sockRedirect")) {
		cmd <- paste(c("proc sockRedirect {sock tosock} {",
			"if {[eof $sock] == 1 || [catch {gets $sock line}]} {",
			"    # end of file or abnormal connection drop",
			"    fileevent $sock readable {}",
			"    close $sock",
			"} else {",
			"    puts $tosock $line",
			"}\n}"),
			collapse = "\n")
		.Tcl(cmd)
	}

	## Create the socket server itself in Tcl (a different one for each port)
	## If we want a secure server, use the tls secured socket instead
	if (isTRUE(secure)) {
		.Tcl(paste("set Rserver_", port, "(main) [tls::socket -server sockAccept_",
			#port, " -require 1 -cafile caPublic.pem -certfile ~/serverR.pem ",
			port, " -certfile Rserver.pem -keyfile Rserver.pem -ssl2 1 -ssl3 1 -tls1 0 -require 0 -request 0 ",
			port, "]", sep =""))
			## For client, use:
			## set chan [tls::socket -cafile caPublic.pem -certfile ~/clientR.pem server.site.net $port]
			## To generate the keys:
			## cd ~
			## Copy /System/Library/OpenSSL/openssl.cnf on ~, and edit
			## openssl genrsa -out serverR.pem 1024   # use -des3 to secure with a password
			## openssl req -new -x509 -key serverR.pem -out clientR.pem -days 365 -config openssl.cnf
			## ... and answer to a couple of questions
	} else {
		.Tcl(paste("set Rserver_", port, "(main) [socket -server sockAccept_",
			port, " ", port, "]", sep =""))
	}

	## Add this port in the variable 'SocketServers' in Sciviews:TempEnv
	socks <- getSocketServers()
	namesocks <- names(socks)
	if (!(portnum %in% socks)) {
		socks <- c(socks, portnum)
		names(socks) <- c(namesocks, server.name)
		socks <- sort(socks)
		assign("SocketServers", socks, envir = TempEnv())
	}
    return(TRUE)  # Humm! Only if it succeeds...
}
