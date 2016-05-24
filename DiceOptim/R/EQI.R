#source("EQI.R")
EQI <- function(x, model, new.noise.var=0, beta=0.9, q.min=NULL, type = "UK", envir=NULL)
{
########## Convert x in proper format(s) ###
d <- length(x)
if (d != model@d){ stop("x does not have the right size") }
newdata.num <- as.numeric(x)
newdata <- data.frame(t(newdata.num))
colnames(newdata) = colnames(model@X)

######### Compute q.min if missing #########
if (is.null(q.min))
{ pred <- predict.km(model, newdata=model@X, type=type, checkNames = FALSE)
  q.min <- min(pred$mean + qnorm(beta)*pred$sd)  }

######### Compute prediction at x #########
predx <- predict.km(model, newdata=newdata, type=type, checkNames = FALSE)
mk.old <- predx$mean
sk.old <- predx$sd

######### Intermediate values ##########
v <- predx$Tinv.c
c <- predx$c

######### mq and sq ##########
if (sk.old < sqrt(model@covariance@sd2)/1e6)
{ mq <- mk.old
  sq <- 0
} else
{ mq <- mk.old + qnorm(beta) * sqrt((new.noise.var * sk.old^2)/(new.noise.var + sk.old^2))
  sq <- sk.old^2/sqrt(new.noise.var + sk.old^2)
}

######### EQI ##########
xcr <- (q.min - mq)/sq
xcr.prob <- pnorm(xcr)
xcr.dens <- dnorm(xcr)

if (!is.null(envir)) {
assign("c", predx$c, envir=envir)
assign("v", v, envir=envir)
assign("xcr.prob", xcr.prob, envir=envir)
assign("xcr.dens", xcr.dens, envir=envir)
assign("mq", mq, envir=envir)
assign("sq", sq, envir=envir)
assign("sk.old", sk.old, envir=envir)
}

if (sk.old < sqrt(model@covariance@sd2)/1e6){ eqi.val <- 0} else 
{ eqi.val <- sq*(xcr*xcr.prob + xcr.dens) }

return(eqi.val)
}
