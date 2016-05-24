parSocket <- function (client, serverport = 8888, clientsocket = client, ...)
{
    ## Set or get parameters for a given socket client
    ## No attempt is made to make sure this client exists
    sc <- paste("SocketClient", client, sep = "_")
    if (!exists(sc, envir = TempEnv(), inherits = FALSE, mode = "environment")) {
        ## Create a new environment with default values
        e <- new.env(parent = TempEnv())
        e$client <- client
		e$clientsocket <- clientsocket
        e$serverport <- serverport
        e$prompt <- ":> "    # Default prompt
        e$continue <- ":+ "  # Default continuation prompt
        e$code <- ""         # Current partial code for multiline mode
        e$last <- ""         # String to add at the end of evaluations
        e$echo <- FALSE      # Don't echo commands to the console
        e$flag <- FALSE      # Do not flag pieces of code (not used yet!)
        e$multiline <- TRUE  # Allow for multiline code
        e$bare <- TRUE       # Always start in "bare" mode
        ## Note: in bare mode, all other parameters are inactive!
        ## Assign it to SciViews:TempEnv
        assign(sc, e, envir = TempEnv())
    } else e <- get(sc, envir = TempEnv(), mode = "environment")
    ## Change or add parameters if they are provided
    ## There is no reason that serverport changes
	## but if a client disconnects and reconnects, the clientsocket may be
	## different! But only change if it is sockXXX
	if (grepl("^sock[0-9]+$", clientsocket)) e$clientsocket <- clientsocket
	args <- list(...)
    if (l <- length(args)) {
        change.par <- function (x, val, env) {
            if (is.null(x)) return(FALSE)  # Do nothing without a valid name
            if (is.null(val)) {
                suppressWarnings(rm(list = x, envir = env))  # Remove it
                return(TRUE)
            }
            env[[x]] <- val  # Add or change this variable in the environment
            return(TRUE)
        }
        n <- names(args)
        res <- rep(TRUE, l)
        for (i in seq_len(l))
			res[i] <- change.par(n[i], args[[i]], e)
        if (any(!res))
			warning("Non named arguments are ignored")
    }
    ## Return e invisibly
    return(invisible(e))
}
