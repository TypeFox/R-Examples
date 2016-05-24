ShortenStrings <- function(strings, max.len=40, end.len=floor(max.len/2), sep='...'){
  # for strings longer than max.len, do the following
  # compute: start.len = max.len - end.len - length(sep)
  # then replace characters between start and end with sep
  stopifnot(max.len >= (1+str_length(sep)) & end.len >= 0)
  if (max.len < (end.len + str_length(sep))){
    # the end length is too long relative to sep and max.len
    # reduce end length to max allowable value
    end.len <- max.len - str_length(sep)
    #stop('end.len must be >= end.len + str_length(sep)')
  }
  
  x <- as.character(strings)
  
  # loop for each string in x
  for (i in 1:length(x)){
    xstring <- x[i]
    if (str_length(xstring) <= max.len){
      # string is already short, skip it
      next
    } else {
      # try removing whitespace from start / end
      xstring <- str_trim(xstring)
    }
    # check length again
    if (str_length(xstring) <= max.len){
      # string is short enough after trim, update and move to next string
      x[i] <- xstring
      next
    }
    
    # string is still too long, replace middle section
    start.len <- max.len - str_length(sep) - end.len
    stopifnot(start.len >= 0)
    # handle the case where the end length is zero
    if (end.len == 0){
      after.sep <- ''
    } else {
      after.sep <- str_sub(xstring, -end.len, -1)
    }
    x[i] <- paste(str_sub(xstring, 1, start.len), sep, after.sep, sep='')
  }
  return(x)
}