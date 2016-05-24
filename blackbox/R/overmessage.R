overmessage <- function(msg, prevmsglength) {
  msglength <- nchar(msg)
  if (prevmsglength>0) {base::message("\r", appendLF=F)}  	##FR: for backslash-b see ?Quotes ...
  base::message(msg, appendLF=F)
  return(msglength)
}

overcat <- function(msg, prevmsglength) {
  msglength <- nchar(msg)
  if (prevmsglength>0) {cat("\r")}  	##FR: for backslash-b see ?Quotes ...
  cat(msg)
  return(msglength)
}
