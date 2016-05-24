# ==> clientid&id="xxxxx"
# <==  code=0
# Sends the server an identification of the client, typically a program name.

clientid <- function(id = paste("seqinr_", packageDescription("seqinr")$Version, sep = ""), socket, verbose = FALSE){
  #
  # Client ID definition : seqinr + package version number
  #  (internal note: log file is: /mnt/users/ADE-User/racnuc/log)
  #
  request <- paste("clientid&id=", id, sep = "")
  if(verbose) cat(paste("clientid(): sending", request, "\n"))
  writeLines( request, socket, sep = "\n")
  rep <- readLines(socket, n = 1)
  if(verbose) cat(paste("... answer from server is:", rep, "\n"))
  res <- parser.socket(rep, verbose = verbose)
  if( res[1] != "0") {
    print(rep)
    stop("I don't know what this error code means for clientid, please contact package maintener.\n")
  }
}

