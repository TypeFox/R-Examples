#' @title define Classes 2
#' @description creates a list with the class interval of a frequency distribution
#' @param dados a vector of numbers
#' @param amplitude integer Class amplitude range
#' @examples
#' dados <- c(1,2,3,4)
#' defineClasses2(dados,2)
#' @export
defineClasses2 = function(dados, amplitude) {
  (pCentro = seq(floor(min(dados)),ceiling(max(dados)), amplitude))

  (pClasse = cbind(as.matrix(pCentro)-(amplitude/2),as.matrix(pCentro+amplitude)-(amplitude/2)))

  if (pClasse[dim(pClasse)[1], 2] == max(dados)) {
    pClasse[dim(pClasse)[1], 2] =  max(dados)+0.00000001
  }

  return(list(centro=pCentro, classe=pClasse))
}
