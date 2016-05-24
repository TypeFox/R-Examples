#source("kriging.quantile.R")

kriging.quantile <- function(x, model, beta=0.1, type = "UK", envir=NULL)
{
  ########## Convert x in proper format(s) ###
  d <- length(x)
  if (d != model@d){ stop("x does not have the right size") }
  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) = colnames(model@X)
  
  # Prediction en newdata en partant de X
  predx <- predict(model, newdata=newdata, type=type, checkNames = FALSE)
  mk <- predx$mean
  sk <- predx$sd  
  qk <- mk + qnorm(beta)*sk
  
  if (!is.null(envir)) {
    assign("mk", mk, envir=envir)
    assign("sk", sk, envir=envir)
    assign("c", predx$c, envir=envir)
    assign("Tinv.c", predx$Tinv.c, envir=envir)
  }
  
  return(res <- qk)
}
