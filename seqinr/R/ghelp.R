# ==>   ghelp&file=xx&item=xx
# <==  nl=xx&...1 or several lines...
# Reads one item of information from specified help file.
# File can be HELP or HELP_WIN, item is the name of the desired help item
# Reply :	nl is 0 if any problem, or announces the number of help lines returned.

ghelp <- function(item = c("GENERAL", "SELECT", "SPECIES", "KEYWORD"), 
                  file = c("HELP", "HELP_WIN"), socket = autosocket(),
                  catresult = TRUE){
  #
  # Default is "HELP" file and "GENERAL":
  #
  item <- item[1]
  file <- file[1]
  if(!(file %in% c("HELP", "HELP_WIN"))) stop("Wrong file agument")
  #
  # Build request:
  #
  request <- paste("ghelp&file=", file, "&item=", item, sep = "")
  writeLines(request, socket, sep = "\n")
  answerFromServer <- readLines(socket)
  #
  # Check that there is an answer from server:
  #
  if(length(answerFromServer) == 0){
    warning("Empty answer from server")
    return(NA)
  }
  #
  # cat result:
  #
  answerFromServer[1] <- unlist(strsplit(answerFromServer[1], split = "&"))[2]
  if(catresult) cat(answerFromServer, sep = "\n")
  invisible(answerFromServer)
}

