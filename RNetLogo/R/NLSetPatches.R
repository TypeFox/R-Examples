NLSetPatches <-
function(patch.var, in.matrix, nl.obj=NULL)
{
  # get internal nl.obj if NULL
  if (is.null(nl.obj))
  {
    nl.obj <- "_nl.intern_"
  }
  # check for unknown nl.obj
  if (!(nl.obj %in% .rnetlogo$objects)) {
    stop(paste('There is no NetLogo reference stored under the name ',nl.obj,".", sep=""))
  }   

  if (.rnetlogo$nl3d == TRUE) {
    stop('NLSetPatches does not work in NetLogo 3D! Use NLSetPatchSet instead.')
  }

  if (!is.matrix(in.matrix))
  {
    stop('First argument must be a matrix!')
  }
  min.pxcor <- NLReport("min-pxcor", nl.obj=nl.obj)
  max.pxcor <- NLReport("max-pxcor", nl.obj=nl.obj)
  min.pycor <- NLReport("min-pycor", nl.obj=nl.obj)
  max.pycor <- NLReport("max-pycor", nl.obj=nl.obj)
  
  xdim <- min.pxcor:max.pxcor
  ydim <- max.pycor:min.pycor
  
  dims <- dim(in.matrix)
  if ((length(xdim) != dims[2]) || (length(ydim) != dims[1]))
  {
    stop(paste('matrix dimensions (',dims[2],dims[1],') do not fit NetLogo World dimensions (',length(xdim),length(ydim),')'))
  }
  
  prev_ <- "(foreach sort patches ["  
  inp_ <- paste(t(in.matrix), collapse=" ")                 
  between_ <- "] "
  ask_ <- "[ask ?1 [ "
  set_ <- paste(" set ", patch.var, " ?2", sep="")
  end_ <- "]])"
  
  merged_ <- paste(prev_, inp_, between_, ask_, set_, end_, sep="")
  NLCommand(merged_, nl.obj=nl.obj)    
}

