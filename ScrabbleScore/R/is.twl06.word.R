data(sysdata,envir=environment())
is.twl06.word <- function(w){
  w <- tolower(w)
  w %in% twl06
}