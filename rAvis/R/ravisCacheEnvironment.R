# cache mecanism for internal use of the package

.ravis_cache <- new.env()

# returns a cached variable if it's intializated. If not, 
# initialize it by calling the function provided
# 
.avisCacheReturnOrSetup<-function(varname, setupfun)
{
  if(!.avisCacheHas(varname))
  {
  	if(is.function(setupfun)){
  		setupcall <- as.call(list(setupfun))
    } else {
      setupcall <- call(setupfun)
    }

    setupvalue <- eval(setupcall)
    .avisCacheSet(varname, setupvalue)

    .avisCacheGet(varname)
  }

  return (.avisCacheGet(varname))
}

.avisCacheHas<-function(varname)
{
	exists(varname, envir=.ravis_cache)
}

.avisCacheGet<-function(varname)
{
	return (get(varname, envir=.ravis_cache))
}

.avisCacheSet<-function(varname, value)
{
	assign(varname, value, envir = .ravis_cache)
}