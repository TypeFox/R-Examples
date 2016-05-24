NLQuit <-
function(nl.obj=NULL, all=FALSE)
{
  if (all) {
    # object names
    #objs <- names(.rnetlogo$objects)
    objs <- .rnetlogo$objects

    # handle gui obj as last, to prevent java.lang.InterruptionException
    guiobj <- NULL    
    if ((!is.null(.rnetlogo$guiobj)) && (.rnetlogo$guiobj %in% objs)) {
      objs <- objs[which(!objs == .rnetlogo$guiobj)]
      guiobj <- .rnetlogo$guiobj
    }
    
    invisible(
      lapply(objs, function(x) {NLQuit(nl.obj=x, all=F)})
    )  
    
    # close gui obj
    if (!is.null(guiobj)) {
      NLQuit(nl.obj=guiobj, all=F)
    }
    
  } else { 
    obj.name <- nl.obj
    if (is.null(obj.name))
    {
      obj.name = "_nl.intern_"
    }
    if (obj.name %in% .rnetlogo$objects) {
	  nl.obj <- get(obj.name, envir=.rnetlogo)
    } else {
      stop(paste('There is no NetLogo reference stored under the name ',obj.name,".", sep=""))
    }    
    
    .jcall(nl.obj, "V", "KillWorkspace")
  
    #print("jinit")
    #.jinit(force.init=TRUE)
    #print("detach RNetLogo")
    #detach("package:RNetLogo")
  
    # see http://osdir.com/ml/lang.r.rosuda.devel/2007-01/msg00052.html
    #print("doneJVM")
    #.Call("doneJVM")
  
    #print("detach rJava")
    #detach("package:rJava")
  
    # free the instance
    # doesn't work for others than .rnetlogo[['nl.intern']]
    nl.obj <- NULL
    #.rnetlogo$objects[obj.name] <- NULL
    .rnetlogo$objects <- .rnetlogo$objects[-which(.rnetlogo$objects %in% obj.name)]
	
    # call the garbage collector
    .jcall('java/lang/System', 'V', 'gc')
    # java error handling
    if (!is.null(e<-.jgetEx()))
    {
      if (.jcheck(silent=TRUE))
      {
        print(e)
        stop()
      }
    } 
    
    # reset working directory after last NetLogo instance was closed
    if (length(.rnetlogo$objects) == 0) {
      setwd(.rnetlogo$savedworkingdir[1])
    }
  }
}

