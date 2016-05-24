residuals.glm1path = function(object, ... )
{
    tol = 1.e-8
    n.rows = length(object$y)
    fits    = object$glm1.best$fitted
    
    if(is.character(object$family))
      familyString = object$family
    else
      familyString = object$family$family
      
    if (is.na(pmatch("Negative",familyString))==FALSE) #if a proper MASS negative binomial family argument was specified:
      familyString = "negative.binomial"
    
    if(familyString=="gaussian") #work out residual variance if gaussian
    {
      df.residual = n.rows - min(object$df[object$lambdas==object$lambda])
      sigma2 = (object$y-fits)^2 / df.residual
    }
    
    pfn = switch(familyString, #set the function to use for cdf calculation
                 "negative.binomial" = "pnbinom",
                 "poisson" = "ppois",
                 "binomial" = "pbinom",
                 "gaussian" = "pnorm",
                 stop("unknown character vector for family argument")  
                )

    params = switch(familyString, #set the parameters for call to pnf
                    "negative.binomial" = list(q=object$y, mu=fits, size=object$glm1$phi),
                    "poisson" = list(q=object$y, lambda=fits),
                    "binomial" = list(q=object$y, size=1, prob=fits ),
                    "gaussian" = list(q=object$y, mu=fits, sd=sqrt(sigma2))
                    )

    param.minus = params
    param.minus$q = params$q - 1.e-6 #for discrete cases, (hopefully!) ignorable for continuous cases

    
    u = runif(n.rows)
      
    qupper = do.call(pfn, params)
    qlower = do.call(pfn, param.minus)
    resids = u * pmax( tol^3, qupper ) + (1-u) * pmin( 1-tol^3, qlower )
      #pmax and pmin used to avoid any values identically 0 or 1

    return( qnorm(resids) )
}


plot.glm1path = function(x, add.smooth=getOption("add.smooth"), ... )
{
  resid      = residuals.glm1path(x)  
  fits       = x$glm1.best$fitted
  linearPred = x$glm1.best$family$linkfun(fits)
  plot(linearPred, resid, xlab="Linear predictor", ylab="Dunn-Smyth residual", type="n",...)
  abline(h=0,col="red",lwd=2)
  points(linearPred,resid, ...)
  if(add.smooth==TRUE)
  {
      smoother=loess.smooth(linearPred,resid)
      lines(list(x=smoother$x,y=smoother$y), col="red")
  }
  points(linearPred, resid, ...)
}

