#' @title calculates percentage
#' @description With this function, you can calculate the ratio of one quantity or magnitude relative to another evaluated in percentage.
#' @param valor number amount you to know the percentage
#' @param observados number relationship to which you want to calculate the percentage, if it is a vector of integers is calculated its average.
#' @details calculaPerc = ((valor)/mean(observados))*100
#' @export
calculaPerc <- function(valor, observados)
{
  ((valor)/mean(observados))*100
}
