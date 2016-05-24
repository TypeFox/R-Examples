###############################################################################
#                                                                             #
#                             parser.socket                                   #
#                                                                             #
#          Utility function to parse answers from ACNUC server.               #
#                                                                             #
###############################################################################

parser.socket <- function(onelinefromserver, verbose = FALSE)
{
  if(verbose) cat(paste("parser.socket received: -->",
                        onelinefromserver,"<--\n", sep = ""))
  if(length(onelinefromserver) == 0){
    if(verbose) cat("character(0) detected returning NULL\n")
    return(NULL)
  }
  if(is.null(onelinefromserver)){
    if(verbose) cat("NULL detected returning NULL\n")
    return(NULL)
  }
  #
  # Answer from server looks like: "code=0&lrank=2&count=150513&type=SQ&locus=F"
  # 
  loc <- gregexpr("=[^=&]*", onelinefromserver)[[1]]
  substring(onelinefromserver, loc + 1, loc + attr(loc, "match.length") - 1)
}
