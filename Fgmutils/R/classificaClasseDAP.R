#' @title classifica Classe DAP
#' @description the center of the class that the DAP belongs.
#' @param dfClassesDAP a frequency distribution with the attributes $classe and $centro
#' @param dap integer Diameter at breast height
#' @examples
#' dados <- c(1,2,3,4)
#' dados = defineClasses2(dados,2)
#' classificaClasseDAP(dados,2)
#' @export
classificaClasseDAP <- function(dfClassesDAP, dap) {
  for (i in 1:length(dfClassesDAP$centro)) {
    classe = dfClassesDAP$classe[i,]
    centro = dfClassesDAP$centro[i]
    linf=classe[1]
    lsup=classe[2]

    if (dap>=linf && dap<lsup) {
      return(centro)
    }
  }
  return(NA)
}
