.metagearPROBLEM <- function(type, aMessage) {
  newMessage <- paste0("metagear ",
                       type,
                       " in ",
                       as.list(sys.call(-1))[[1]],
                       "(): ",
                       aMessage,
                       ".")
  if(type == "error") stop(newMessage, call. = FALSE)
  message(newMessage)
}

renameFile <- function(aFileName) {
  fileNumber <- gsub("[^0-9]", "", aFileName)
  aFileName <- ifelse(fileNumber == "", 
                      gsub("\\.", "1.", aFileName),
                      gsub("\\d+", paste0(as.integer(fileNumber) + 1), aFileName))
  return(aFileName)  
}