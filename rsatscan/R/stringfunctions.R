stripslash = function (string) {
  laststr = tail(unlist(strsplit(string,"")),1)
  return( ifelse(laststr == "/", substr(string, 1, nchar(string)-1) , string)) 
}

