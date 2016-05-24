#' Calculate VIP scores for PLS regression
#' 
#' This function calculates the Variable Importance in the Projection statistic
#' for the Partial Least Squares regression. It is used in the PLS function.
#' Executing it in isolation will probably not be useful to most users.
#' 
#' This is required to produce the VIP scores for the PLS procedure.
#' 
#' @param object an mvr object, as produced by the pls procedure or a range of
#' other functions
#' @return data frame with as many columns as independent variables are input
#' into the PLS analysis. The number of columns corresponds to the number of
#' latent components selected for the analysis. Values in the data frame are
#' the VIP values corresponding to each variable for the respective component.
#' @author Eike Luedeling, but the function was mainly copied from
#' http://mevik.net/work/software/pls.html; the reference given there is listed
#' below
#' @references the function is mostly identical to the one provided on
#' http://mevik.net/work/software/pls.html.
#' 
#' Here is the reference given there:
#' 
#' Chong, Il-Gyo & Jun, Chi-Hyuck, 2005, Performance of some variable selection
#' methods when multicollinearity is present, Chemometrics and Intelligent
#' Laboratory Systems 78, 103-112
#' 
#' This reference refers to the chillR package:
#' 
#' Luedeling E, Kunz A and Blanke M, 2013. Identification of chilling and heat
#' requirements of cherry trees - a statistical approach. International Journal
#' of Biometeorology 57,679-689.
#' @keywords utility
#' @examples
#' 
#' PLS_results<-PLS_pheno(
#'   weather_data=KA_weather,
#'   split_month=6,   #last month in same year
#'   bio_data=KA_bloom,return.all=TRUE)
#' 
#' #return.all makes the function return the whole PLS object - needed for next line to work
#'   
#' VIP(PLS_results$PLS_output)
#'   
#'  
#' @export VIP
VIP <-
function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}
