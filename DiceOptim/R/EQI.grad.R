EQI.grad <- function(x, model, new.noise.var=0, beta=0.9, q.min=NULL, type = "UK", envir=NULL){

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

######### Intermediate values ##########
T <- model@T
z <- model@z
u <- model@M

if (is.null(envir))
{
  ######### Compute prediction at x #########
  predx <- predict(model, newdata=newdata, type=type, checkNames = FALSE)
  mk.old <- predx$mean
  sk.old <- predx$sd
  v <- predx$Tinv.c
  c <- predx$c
  
  ######### mq and sq ##########
  if (sk.old < 1e-16)
  { mq <- mk.old
    sq <- 0
  } else
  { mq <- mk.old + qnorm(beta) * sqrt((new.noise.var * sk.old^2)/(new.noise.var + sk.old^2))
    sq <- sk.old^2/sqrt(new.noise.var + sk.old^2)
  } 
  
  xcr <- (q.min - mq)/sq
  xcr.prob <- as.numeric(pnorm(xcr))
  xcr.dens <- as.numeric(dnorm(xcr))
  
} else
{ toget <- matrix(c("xcr.prob", "xcr.dens", "c", "v", "mq", "sq"),1,12)
  apply(toget, 2, get, envir=envir)
  xcr.prob <- envir$xcr.prob
  xcr.dens <- envir$xcr.dens
  c        <- envir$c
  v        <- envir$v
  mq       <- envir$mq
  sq       <- envir$sq
  sk.old   <- envir$sk.old
}

if (sk.old < sqrt(model@covariance@sd2)/1e6)
{ eqi.grad.val <- rep(0,d)
} else
{
  ######### Intermediate values ##########
  F.newdata <- model.matrix(model@trend.formula, data=newdata)
  dc <- covVector.dx(x=newdata.num, X=model@X, object=model@covariance, c=c)
  f.deltax <- trend.deltax(x=newdata.num, model=model)
  W <- backsolve(t(T), dc, upper.tri=FALSE)
  
  ######### Gradient of mk and sk2 of the old model ##########
  mk.old.grad <- t(W)%*%z + t(model@trend.coef%*%f.deltax)
  if (type=="UK")
  { tuuinv <- solve(t(u)%*%u)
    sk2.old.grad <-  t( -2*t(v)%*%W + 2*(F.newdata - t(v)%*%u )%*% tuuinv %*% (f.deltax - t(t(W)%*%u) ) )
  } else
  { sk2.old.grad <-  t( -2*t(v)%*%W) }
  
  ######### Gradients of mq and sq ##########
  mq.grad  <- mk.old.grad + as.numeric(qnorm(beta)*(new.noise.var)^(3/2)/(2*sk.old*(sk.old^2 + new.noise.var)^(3/2)))*sk2.old.grad
  sq2.grad <- as.numeric((sk.old^2)*(2*new.noise.var + sk.old^2)/(new.noise.var + sk.old^2)^2)*sk2.old.grad
  sq.grad <- sq2.grad/as.numeric(2*sq)
  
  ######### Gradient of EQI ##########
  eqi.grad.val <- - xcr.prob * mq.grad + xcr.dens * sq.grad
}
return(eqi.grad.val)
}
