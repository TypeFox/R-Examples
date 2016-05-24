# ==>  knowndbs{&tag=xx}
# <== nl=.. \n
# dbname | on/off | db description \n      nl such lines
# Returns, for each database known by the server, its name (a valid value for the db= argument 
# of the acnucopen command), availability (off means temporarily unavailable), and description.
# When the optional tag= argument is used, only databases tagged with the given string are listed;
# without this argument, only untagged databases are listed.
# The tag argument thus allows to identify series of special purpose (tagged) databases, 
# in addition to default (untagged) ones. The full list of untagged and tagged databases is here.

knowndbs <- function(tag = c(NA, "TP", "TEST", "DEV"), socket = autosocket()){
  #
  # Use default tag:
  #
  tag <- tag[1]
  #
  # Build request:
  #
  if( !is.na(tag) ){
    askforbank <- paste("knowndbs&tag=", tag, sep = "")
  } else {
    askforbank <- "knowndbs"
  }
  writeLines(askforbank, socket, sep = "\n")
  rep <- readLines(socket, n = 1)
  nbank <- as.numeric(parser.socket(rep))  
  #
  # Read bank infos from server:
  #
  res <- readLines(socket, n = nbank)
  #
  # Build result:
  #
  resdf <- as.data.frame(list(bank = I(rep("NAbank", nbank)), 
                  status = I(rep("NAstatus", nbank)), 
                  info = I(rep("NAinfo", nbank))))
  for(i in seq_len(nbank))
    resdf[i, ] <- unlist(strsplit(res[i], split = "\\|"))[1:3]
  for(i in seq_len(nbank))
    for(j in seq_len(3))
      resdf[i, j] <- trimSpace(resdf[i, j])   
  #
  # Return result:
  #
  return(resdf)
}

kdb <- knowndbs
