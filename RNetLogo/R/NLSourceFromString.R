NLSourceFromString <-
function(..., append.model=TRUE, nl.obj=NULL)
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
  
  if(.rnetlogo$startedGUI)
  {
    # preprocessing: evaluate the commands
    model.source <- lapply(list(...), function(x) {eval.commandobject(x)})
    # put all commands together to one string to be evaluated by NetLogo
    model.source <- paste(c(...), collapse="\n")
  
    .jcall(nl.obj, "V", "sourceFromString", .jnew("java/lang/String", model.source), .jnew("java/lang/Boolean",append.model))
    # java error handling
    if (!is.null(e<-.jgetEx()))
    {
      if (.jcheck(silent=TRUE))
      {
        print(e)
        stop()
      }
    }
  }
  else
  {
    stop('NLSourceFromString is only available if NetLogo was started with GUI.')
  }
}

