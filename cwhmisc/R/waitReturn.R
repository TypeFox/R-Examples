waitReturn <- function(q="",ask=TRUE) {
  if (ask & interactive() & sink.number()==0) readline(paste(q,"\nType  <Return>\t to continue : "))
  invisible()
}
