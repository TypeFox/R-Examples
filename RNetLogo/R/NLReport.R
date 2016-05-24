NLReport <-
function(reporter, nl.obj=NULL)
{
  # get internal nl.obj if NULL
  if (is.null(nl.obj))
  {
    nl.obj <- "_nl.intern_"
  }
  # get NetLogo reference
  if (nl.obj %in% .rnetlogo$objects) {
	nl.obj <- get(nl.obj, envir=.rnetlogo)	
  } else {
    stop(paste('There is no NetLogo reference stored under the name ',nl.obj,".", sep=""))
  } 
    
  if (length(reporter) == 1) 
  {
    resobj <- .jcall(nl.obj, "Ljava/lang/Object;", "report", .jnew("java/lang/String", reporter))       
    resobj <- eval.reportobject(resobj)
  }
  else {
    resobj <- .jcall(nl.obj, "[Ljava/lang/Object;", "report", .jarray(reporter))
    resobj <- lapply(resobj, function(x) {eval.reportobject(x)})
  } 
  # java error handling
  if (!is.null(e<-.jgetEx()))
  {
    if (.jcheck(silent=TRUE))
    {
      print(e)
      stop()
    }
  } 
  return (resobj)
}

