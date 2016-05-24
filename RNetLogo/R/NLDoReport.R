NLDoReport <-
function(iterations, command, reporter, as.data.frame=FALSE, df.col.names=NULL, nl.obj=NULL)
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
    resobj <- .jcall(nl.obj, "[Ljava/lang/Object;", "doReport", 
                    .jnew("java/lang/String", command), 
                    .jnew("java/lang/String", reporter), 
                    .jnew("java/lang/Integer", as.integer(iterations))
                    )       
    resobj <- lapply(resobj, function(x) {eval.reportobject(x)})
  } 
  else 
  {
    resobj <- .jcall(nl.obj, "[Ljava/lang/Object;", "doReport", 
                    .jnew("java/lang/String", command), 
                    .jarray(reporter), 
                    .jnew("java/lang/Integer", as.integer(iterations))
                    )       
    resobj <- lapply(resobj, function(x) {.jevalArray(x)})
    resobj <- lapply(resobj, function(x) {
                        lapply(x, function(z) {eval.reportobject(z) } )
                      }
                    )
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
  
  # transform result to data.frame
  if (as.data.frame == TRUE)
  {
    resobj <- data.frame(do.call("rbind",resobj))
    length.of.col <- lapply(resobj, function(x) {length(x[[1]])})
    resobj[c(which(length.of.col==1))] <- as.data.frame(lapply(resobj[c(which(length.of.col==1))], function(x) { unlist(x) })) 
    if (length(df.col.names) > 0)
    {
      names(resobj) <- df.col.names
    }
  }
  return (resobj)
}

