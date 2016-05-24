##' Convert a character pattern into NA in character and vice versa.
##'
##' Convert a character pattern into NA in character and vice versa.
##' @param data a data.frame to summarize.
##' @param na.character a character vector specifying missing character.
##' @param reverse a logical value specifying reverse replacement that NA is replaced with the first element of na.character.
##' @return a data.frame with NA replacement.
##' @seealso
##' \code{\link{countNA}}
##' @examples
##' dat <- data.frame(a=c("","2","3"),b=c("4", NA, "."), stringsAsFactors=FALSE)
##' dat2 <- convertNA(dat)
##' dat3 <- convertNA(dat2, na.character=".", reverse=TRUE)
##' dat
##' dat2
##' dat3
##' @export
convertNA <- function(data=NULL, na.character=c("", "."), reverse=FALSE){

  if (is.null(data) || !is.data.frame(data)) 
    stop("data is not specified or not data.frame")
  
  n <- nrow(data)
  
  for(i in 1:ncol(data)){
    if(!is.character(data[,i]))
      next
    for(j in 1:n){
      if(reverse){
        if(is.na(data[j,i]))
           data[j,i] <- na.character[1]
         } else {
        if(data[j,i] %in% na.character)
          data[j,i] <- NA
      }
    }
  }
  
  return(data)
}
