#####################################################################################
# returns the residuals of a manyglm #
#####################################################################################

residuals.manyglm<- function(object, ...)
{
  tol=1.e-8
  n.rows = NROW(object$y)
  n.vars = NCOL(object$y)
  params = list()
  
  if(object$family=="negative.binomial")
  {
    pfn = "pnbinom"
    for(i.var in 1:n.vars)
      params[[i.var]]=list(q=object$y[,i.var],mu=object$fitted[,i.var],size=object$theta[i.var])
  } 
  else if(object$family=="poisson")
  {
      pfn = "ppois"
      for(i.var in 1:n.vars)
        params[[i.var]]=list(q=object$y[,i.var],lambda=object$fitted[,i.var])
  }    
  else if (substr(object$family,1,3) == "bin" || "clo")
  {
    pfn = "pbinom"
    for(i.var in 1:n.vars)
      params[[i.var]]=list(q=object$y[,i.var], size=1, prob=object$fitted[,i.var] )
  }
  else if(object$family=="gaussian")
  {
    pfn = "pnorm"
    df.residual = n.rows - dim(coef(object))[1]
    sigma2 = apply( (object$y-object$fitted)^2, 2, sum ) / df.residual
    for(i.var in 1:n.vars)
      params[[i.var]] = list(q=object$y[,i.var], mu=object$fitted[,i.var], sd=sqrt(sigma2[i.var]))
  } 
  else stop (paste("'family'", object$family, "not recognized"))
  
  resids=matrix(NA,n.rows,n.vars)
  dimnames(resids)[[1]] = rownames(object$y)
  dimnames(resids)[[2]] = names(object$y)
  for(i.var in 1:n.vars)
  {
    param.minus = params[[i.var]]
    param.minus$q = params[[i.var]]$q - 1.e-6
    u = runif(n.rows)

    qupper = do.call(pfn, params[[i.var]])
    qlower = do.call(pfn, param.minus)
    resids[,i.var] = u * pmax( tol^3, qupper ) + (1-u) * pmin( 1-tol^3, qlower )
    #pmax and pmin used to avoid any values identically 0 or 1
  } 
  return( qnorm(resids) )
}