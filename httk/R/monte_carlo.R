# This function performs Monte Carlo to assess uncertainty and variability for pharmacokinetic models.
monte_carlo <- function(params,which.quantile=0.95,cv.params=NULL,censored.params=NULL,samples=1000,name.model='calc_analytic_css',output.col.model=NA,return.samples=F,...)
{
# A matrix where each column corresponds to a parameter that is varied
# and each row is a different draw of parameter values:  
  MC.matrix <- matrix(NA,nrow=samples,ncol=(length(cv.params)+length(censored.params)))
  colnames(MC.matrix) <- c(names(cv.params),names(censored.params))
# Number of different results to be obtained for the different paramters:
  sample.vec <- rep(NA,samples)
# Any parameter given in censored params is sampled from a normal distribution truncated at zero:
  for (this.param in names(cv.params))
  {
    if (!(this.param %in% names(params))) stop(paste("Cannot find cv.params parameter",this.param,"in parameter list."))
    if (params[[this.param]]>0) MC.matrix[,this.param] <- rtnorm(samples,mean=params[[this.param]],sd=params[[this.param]]*cv.params[[this.param]],lower=0)
    else 
    {
      MC.matrix[,this.param] <- 0  
      warning(paste(this.param,"has mean of zero, yielding SD of zero for fixed cv.  Parameter value fixed at zero."))
    }
  }
# Any parameter given in censored params is sampled from a censored distribution:
  for (this.param in names(censored.params))
  {
    if (!(this.param %in% names(params))) stop(paste("Cannot find censored.params parameter",this.param,"in parameter list."))
    if (!("cv" %in% names(censored.params[[this.param]]))) stop(paste("cv (coefficient of variation) must be specified for parameter",this.param))
    if (!("lod" %in% names(censored.params[[this.param]]))) stop(paste("lod (limit of detection) must be specified for parameter",this.param))
    
    MC.matrix[,this.param] <- r.left.censored.norm(samples,mean=params[[this.param]],sd=params[[this.param]]*censored.params[[this.param]]$cv,lod=censored.params[[this.param]]$lod)
  }
# these.params starts with the default paramter values from the argument 
# provided to the function
  these.params <- params
  for (this.sample in 1:samples)
  {
# Those paramters that have been varied in the Monte Carlo simulation are
# overwritten:  
    these.params[colnames(MC.matrix)] <- MC.matrix[this.sample,]
# Arguments for the do.call statement include "parameters" and anything extra
# given in ...
    these.args <- c(list(parameters = these.params),list(...))
# If model.output.col is set we expect the call to the model to return a list
# with multiple values and we only keep the one in model.output.col
    if (is.na(output.col.model)) 
    {
      out <- do.call(name.model,these.args)
      if (length(out) > 1) stop(paste("Must specific output.col.model because model returned the following arguments:",names(out),collapse=" "))
      sample.vec[this.sample] <- do.call(name.model,these.args)
    } else sample.vec[this.sample] <- do.call(name.model,these.args)[[output.col.model]]
  }
   if(return.samples) out <- sample.vec
   else out <- quantile(sample.vec,which.quantile)

  return(out)
}


