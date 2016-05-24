#' @export lrtest.pers
#' @title Likelihood Ratio Test for Object of class "pers"
#' @description Function to perform a likelihood ratio test for the estimated model 'against' the saturated model for object of class\code{"pers"}. 
#' @param object an object of class\code{"pers"} - see function \code{\link{pers}}.
#' @param ... not used jet.
lrtest.pers<-function(object, ...){
  log_estmod <- sum(unique(object$pers$WLL)) # direktes auslesen der LL aus object
  #log_estmod <- (logLik.pers(object))[[1]]
  df_estmod <- attr(x=logLik.pers(object),which = "df")
  log_satmod <- (logLik.pers(object, sat = TRUE))[[1]]
  df_satmod <- attr(x=logLik.pers(object, sat = TRUE),which = "df")
  chi_test <- max(0,(-2*(log_estmod - log_satmod))) ### 
  df_test <- (df_satmod-df_estmod)
  p_test <- 1-pchisq(q=chi_test, df=df_test)
  result <- list("Chi^2"=chi_test, df = df_test, p = p_test )
  return(result)
}

