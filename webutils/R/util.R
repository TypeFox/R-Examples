# Override default for call. argument
stop <- function(..., call. = FALSE){
  base::stop(..., call. = FALSE)
}

# Strip trailing whitespace
trail <- function(str){
  str <- sub("\\s+$", "", str, perl = TRUE);
  sub("^\\s+", "", str, perl = TRUE);
}

# Remove double quotes if any
unquote <- function(str){
  len <- nchar(str)
  if(substr(str, 1, 1) == '"' && substr(str, len, len) == '"'){
    return(substr(str, 2, len-1));
  } else {
    return(str)
  }
}
