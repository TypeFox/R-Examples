# source("AEI.R")
AEI <- function(x, model, new.noise.var=0, y.min=NULL, type = "UK", envir=NULL)
{
  ########## Convert x in proper format(s) ###
  d <- length(x)
  if (d != model@d){ stop("x does not have the right size") }
  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) = colnames(model@X)

  # Compute y.min if missing
  if (is.null(y.min))
  {
    pred <- predict(object=model, newdata=model@X, type=type, checkNames = FALSE)
    mk <- pred$mean
    sk <- pred$sd
    qk <- mk + qnorm(0.75)*sk
    y.min <- mk[which.min(qk)]
  }
  
  # Prediction en newdata en partant de X
  pred <- predict(object=model, newdata=newdata, type=type, checkNames = FALSE) 
  mk <- pred$mean
  sk <- pred$sd  
  
  xcr <- (y.min - mk)/sk 
  xcr.prob <- pnorm(xcr)
  xcr.dens <- dnorm(xcr)
  
  if (!is.null(envir)) 
  {	assign("xcr", xcr, envir=envir)
    assign("xcr.prob", xcr.prob, envir=envir)
    assign("xcr.dens", xcr.dens, envir=envir)
    assign("c", pred$c, envir=envir)
    assign("Tinv.c", pred$Tinv.c, envir=envir)
    assign("mk", mk, envir=envir)
    assign("sk", sk, envir=envir)
  }
  if (sk < sqrt(model@covariance@sd2)/1e6) { aei.val <- 0 }  else 
  { aei.val <- ((y.min - mk) * xcr.prob + sk * xcr.dens) * (1- sqrt(new.noise.var)/sqrt(new.noise.var + sk^2))}

  return(aei.val)
}
