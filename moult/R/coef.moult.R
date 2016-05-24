coef.moult <- function (object, model = c("full", "duration", "mean", "sd"), ...) 
{ model <- match.arg(model)
  rval <- object$coefficients
  
  rval <- switch(model, 
                 full = structure(c(rval$duration, rval$mean, rval$sd), 
                   .Names = c(paste("duration", names(rval$duration), sep = "_"), 
                     paste("mean", names(rval$mean), sep = "_"),
                     paste("sd", names(rval$sd), sep = "_"))),
                 duration = rval$duration,
                 mean = rval$mean,
                 sd = rval$sd)
  rval
}


