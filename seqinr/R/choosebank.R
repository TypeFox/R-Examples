###########################################################################
# 
#                                choosebank
#
# To select an ACNUC database or to get the list of available databases
# from an ACNUC server.
# 
###########################################################################

.seqinrEnv <- new.env()
choosebank <- function(bank = NA,
                       host = "pbil.univ-lyon1.fr",
                       port = 5558,
                       server = FALSE,
                       blocking = TRUE,
                       open = "a+",
                       encoding = "",
                       verbose = FALSE,
                       timeout = 5,
                       infobank = FALSE,
                       tagbank = NA){
  #
  # Print parameter values if verbose mode is on:
  #
  if(verbose){ 
    cat("Verbose mode is on, parameter values are:\n")
    cat(paste("  bank = ", deparse(substitute(bank)), "\n"))
    cat(paste("  host = ", deparse(substitute(host)), "\n"))
    cat(paste("  port = ", deparse(substitute(port)), "\n"))
    cat(paste("  timeout = ", deparse(substitute(timeout)), "seconds \n"))
    cat(paste("  infobank = ", deparse(substitute(infobank)), "\n"))
    cat(paste("  tagbank = ", deparse(substitute(tagbank)), "\n"))
  }
  
  #
  # Check parameter values (to be completed):
  #
  if( !is.na(tagbank) ){
    if(verbose) cat("I'm checking the tagbank parameter value...\n")
    if( !(tagbank %in% c("TEST", "TP", "DEV")) ){
      if(verbose) cat("... and I was able to detect an error.\n")
      stop("non allowed value for tagbank parameter.\n")
    } else {
      if(verbose) cat("... and everything is OK up to now.\n")
    }  
  }
  
  #
  # Check that sockets are available:
  #
  if(verbose) cat("I'm ckecking that sockets are available on this build of R...\n")
  if( !capabilities("sockets") ){
    stop("Sockets are not available on this build of R.")
   } else {
    if(verbose) cat("... yes, sockets are available on this build of R.\n")
  }
  
  # 
  # Try to open socket connection:
  #
  if(verbose) cat("I'm trying to open the socket connection...\n")
  oldtimeout <- getOption("timeout")
  options(timeout = timeout)
  socket <- try( socketConnection( host = host, port = port, server = server,
                                  blocking = blocking, open = open, encoding = encoding))
  options(timeout = oldtimeout)
  if(inherits(socket, "try-error")) {
    errmess <- paste("I wasn't able to open the socket connection:\n",
                     "  o Check that your are connected to the internet.\n",
                     "  o Check that port", port, "is not closed by a firewall.\n",
                     "  o Try to increase timeout value (current is", timeout, "seconds).\n")
    stop(errmess)
  } else {
    if(verbose) cat("... yes, I was able to open the socket connection.\n")
  }

  #
  # Read the answer from server:
  #
  if(verbose) cat("I'm trying to read answer from server...\n")
  rep1 <- readLines(socket, n = 1)
  if(verbose) cat(paste("... answer from server is:", rep1, "\n"))

  #
  # Send client ID to server:
  #
  clientid(socket = socket, verbose = verbose)
           
  ###############################################################################
  #
  # If no bank name is given, return the list of available banks from server:
  #
  ###############################################################################
  resdf <- kdb(tag = tagbank, socket = socket)
  nbank <- nrow(resdf)
  
  if( is.na(bank) ){
    close(socket) # No more needed
    if(verbose) cat("No bank argument was given...\n")
    if( !infobank ){
      if(verbose) cat("infobank parameter is FALSE, I'm just returning bank names\n")
      return(resdf$bank)
    } else {
      if(verbose) cat("infobank parameter is TRUE, I'm returning all bank infos\n")
      return(resdf) 
      }
  } else {

    ###############################################################################
    #
    # If a bank name is given, try to open it from server:
    #
    ###############################################################################
    
    # 
    # Try to open bank from server:
    #
    if(verbose) cat("I'm trying to open the bank from server...\n")
    resacnucopen <- acnucopen(bank, socket)
    if(verbose) cat("... and everything is OK up to now.\n")

    #
    # Try to get informations from HELP file: 
    #
    if(verbose) cat("I'm trying to get information on the bank...\n")
    bankhelp <- ghelp(item = "CONT", file = "HELP", socket = socket, catresult = FALSE)
    bankrel <- bankhelp[2]
    if(verbose) cat("... and everything is OK up to now.\n")

    #
    # Try to get status info:
    #
    status <- "unknown"
    for(i in seq_len(nbank)){
     if (resdf[i,1] == bank) status <- resdf[i,2]
    }
    
    #
    # Build result and assign it in the global environment:
    #
    res <- list(socket = socket,
     bankname = bank,
     banktype = resacnucopen$type,
     totseqs = resacnucopen$totseqs,
     totspecs = resacnucopen$totspecs,
     totkeys = resacnucopen$totkeys,
     release = bankrel,
     status = status,
     details = bankhelp)
     assign("banknameSocket", res, .seqinrEnv)
    invisible(res)
  }
} 
