##' @title Field Converts To Character
##' @description converts a column of a dataframe to String
##' @param nomeCampo the column name you want to convert
##' @param base the column having dataFrame, that you want to convert to String
##' @return base dataFrame with a column converted to String
##' @examples
##' measurement_date <- c(02/2009,02/2010,02/2011,02/2011)
##' plot <- c(1,2,3,4)
##' test <- data.frame(measurement_date,plot)
##' converteCampoParaCharacter("measurement_date",test)
##' @export
converteCampoParaCharacter <- function (nomeCampo, base){
  eval(parse(text=paste0("base$",nomeCampo," = as.character(base$", nomeCampo, ")")))
  return (base);
}
