###################################################################################################
#                                                                                                 #
#                                         getType                                                 #
#                                                                                                 #
#                                                                                                 #
#                      To get available subsequence types in an opened ACNUC database             #
#                                                                                                 #
###################################################################################################

getType <- function(socket = autosocket()){
  #
  # Get the number of records in SMJYT index file:
  #
  nl <- readfirstrec(type = "SMJ")
  
  #
  # Read the SMJYT index file:
  #
  
  smj <- readsmj(libel.add = TRUE, sname.add = TRUE, nl = nl)
  
  #
  # Return available type:
  #
  ntype <- sum(smj$nature == "type", na.rm = TRUE)
  if( is.na(ntype) || ntype == 0 ){
    return(NA)
  } else {
    return( smj[!is.na(smj$nature) & smj$nature == "type", c("sname","libel")] ) 
  }
}
