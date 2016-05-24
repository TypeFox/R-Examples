xU <- function(x) {

  x <- paste(toupper(substr(x,1,1)), substr(x,2,nchar(x)), sep="")

  return(x)

}
