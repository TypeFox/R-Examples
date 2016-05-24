#' @title return value
#' @description this feature is designed to fix variables that its content was a command
#' @param valor any variable
#' @return the variable converted to its value
#' @examples
#' a = 5
#' retornaValor(a)
#' @export
retornaValor <- function(valor){
  if(typeof(valor) == "character")
  {

    return (as.character(valor));
  }
  return (as.numeric(valor));

}
