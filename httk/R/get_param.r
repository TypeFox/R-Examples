get_param <- function(param.name,Params,calling.func,default=NULL)
{
  if (is.null(Params[[param.name]]))
  {
    if (is.null(default)) stop(paste("Parameter",param.name,"needed in",calling.func))
    else return(default)
  } else return(Params[[param.name]])
}