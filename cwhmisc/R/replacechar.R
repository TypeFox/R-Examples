replacechar <- function(str, char = "_", newchar = ".")
  {  ## tjoelker@redwood.rt.cs.boeing.com (Rod Tjoelker 865-3197)
  under <- grep(char, str)
  for(i in under) {
    nc <- nchar(str[i])
    ch <- substring(str[i], 1:nc, 1:nc)
    ch <- ifelse(ch == char, newchar, ch)
    str[i] <- paste(ch, collapse = "")
  }
  return(str)
}
