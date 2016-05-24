#' @title mspr
#' @description average square of the prediction errors .
#' @param observados vector of values observed.
#' @param estimados vector of regression model data.
#' @param nValidacao number of cases in the validation data set.
#' @references JESUS, S. C.; MIURA, A. K. Analise de regressao linear multipla para estimativa do indice de vegetacao melhorado (EVI) a partir das bandas 3 4 e 5 do sensor TM/Landsat 5. In: SIMPOSIO BRASILEIRO DE SENSORIAMENTO REMOTO, 14. (SBSR), 2009, Natal. Anais... Sao Jose dos Campos: INPE, 2009. p. 1103-1110. DVD, On-line. ISBN 978-85-17-00044-7. (INPE-15901-PRE/10511)
#' @export
mspr <- function(observados, estimados, nValidacao)
{
  sum((observados-estimados))^2/nValidacao
}
