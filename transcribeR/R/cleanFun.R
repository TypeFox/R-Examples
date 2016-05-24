cleanFun <- function(string){
    # Internal function to strip tags from string.
    #    
    # Args:
    #   string: a character string with tags
    #    
    # Returns:
    #   original string with tags removed
  return(gsub("<.*?>", "", string))
}
