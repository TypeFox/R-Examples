#' remove non-ASCII characters
#' @param string any character string
#' @return the string after dropping all html tags to spaces
#' @keywords internal
drop_nonascii <- function(string){
#  string <- gsub("\302\260", " degrees ", string)
#  string <- gsub("\342\200\223", " - ", string)  # <c2><96>
#  Other violating codes <c2><b0>, <c3><a9>, <c3><ad>, <c2><bd>, <c3><b3>
#  Faster solution to drop the offending codes:
  if(is.null(string))
    string <- NULL
  else{
   Encoding(string) <- "UTF-8"
#  string <- iconv(string, "UTF-8", "ASCII", "byte") # return hex code
   string <- iconv(string, "UTF-8", "ASCII", "") # just drop those chars
  }
  string
}

#' clean the fish.base data into pure ASCII
#' @param a list item with fishbase data
#' @return the item scrubbed of non-ASCII characters
#' @keywords internal
clean_data <- function(metadata) lapply(metadata, function(x) lapply(x, drop_nonascii)) 

