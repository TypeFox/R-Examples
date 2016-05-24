#' @title Check de type of Column
#' @description this function returns the type of a column of a dataFrame, if it is numeric or character.
#' @param coluna column of dataframe
#' @examples
#' ID_REGIAO <- c(1,2,3,4)
#' CD_PLANTIO <- c("ACD","CDB","CDC","CDD")
#' test <- data.frame(ID_REGIAO,CD_PLANTIO)
#' verificaTipoColuna(test$ID_REGIAO)
#' @export
verificaTipoColuna <- function(coluna){
  if(typeof(coluna) == "character")
  {

    return ("as.character()");
  }
  return ("as.numeric()");

}
