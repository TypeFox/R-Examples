NLSetPatchSet <-
function(patch.var, input, nl.obj=NULL)
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
  
  if (!is.data.frame(input)) {
    stop("Input has to be a data.frame.")
  }
  
  if (length(patch.var) != 1) {
    stop("You have to submit one patch variable name with argument patch.var")
  }
    
  
  if (.rnetlogo$nl3d == TRUE) {
      # for NetLogo 3D:
    
      if (ncol(input) != 4) {
        stop("Input must have four columns: pxcor, pycor, pzcor, <patch.var>")
      }
      
      start_ <- "(foreach "
      xcords_ <- paste("[", paste(input[['pxcor']], sep=" ", collapse=" "), "]")
      ycords_ <- paste("[", paste(input[['pycor']], sep=" ", collapse=" "), "]")
      zcords_ <- paste("[", paste(input[['pzcor']], sep=" ", collapse=" "), "]")
      var_ <- paste("[", paste(input[[patch.var]], sep=" ", collapse=" "), "]")
      ask_ <- paste("[ ask patch ?1 ?2 ?3 [ set ", patch.var," ?4", sep="") 
      end_ <- " ]])"
      merged_ <- paste(start_, xcords_, ycords_, zcords_, var_, ask_, end_, sep="") 
      NLCommand(merged_, nl.obj=nl.obj)
  }
  else {
    # for conventional 2D NetLogo:      

    if (ncol(input) != 3) {
      stop("Input must have three columns: pxcor, pycor, <patch.var>")
    }
    start_ <- "(foreach "
    xcords_ <- paste("[", paste(input[['pxcor']], sep=" ", collapse=" "), "]")
    ycords_ <- paste("[", paste(input[['pycor']], sep=" ", collapse=" "), "]")
    var_ <- paste("[", paste(input[[patch.var]], sep=" ", collapse=" "), "]")
    ask_ <- paste("[ ask patch ?1 ?2 [ set ", patch.var," ?3", sep="") 
    end_ <- " ]])"
    merged_ <- paste(start_, xcords_, ycords_, var_, ask_, end_, sep="")  
    NLCommand(merged_, nl.obj=nl.obj)
  }
    
}
