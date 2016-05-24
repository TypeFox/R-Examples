NLStart <-
function(nl.path, gui=TRUE, nl.obj=NULL, is3d=FALSE)
  {
  if (is.null(nl.obj))
  {
    nl.obj = "_nl.intern_"
  }
    
  # get RNetLogo jar
  localjar <- 'RNetLogo.jar'
  
  # check for 3D/2D session
  if (.rnetlogo$nl3d == -1)
  {
      .rnetlogo$nl3d <- is3d
  }    
  else
  {
    if (.rnetlogo$nl3d != is3d)
      stop("You can't use 2D and 3D NetLogo in one R session.")
  }
  
  if (!is.null(nl.obj) && (!is.character(nl.obj))) stop("Argument nl.obj has to be a character.")
  if (nl.obj %in% names(.rnetlogo$objects)) stop("Name of object (nl.obj) to store the NetLogo instance
      is already in use. Please quit the used object first or choose a different name.")
      
  # turn off awt in headless mode
  #if (!gui) Sys.setenv(NOAWT=1)
  #else Sys.unsetenv("NOAWT")
  
  
  # NetLogo version check
  if (.rnetlogo$nlversion == 0)
  {
      .rnetlogo$nlversion <- 5
      # load the RNetLogo jar
      .jpackage(
          .rnetlogo$pkgname,
          lib.loc = .rnetlogo$libname,
          jars=localjar
      ) 
  }       
          
  if ((gui) && (.rnetlogo$startedGUI)) stop("RNetLogo was already started with GUI. 
  It isn't possible to start it again in this R session.")
  
  # store working directory, if this is the first NetLogo instance/there are no others
  if (length(.rnetlogo$objects) == 0) {
    .rnetlogo$savedworkingdir <- Prepro(nl.path)
  }    
  
  # use the connection for NetLogo version
  nlo <- tryCatch(
			.jnew("nlcon/NLink",.jnew("java/lang/Boolean",gui),.jnew("java/lang/Boolean",is3d),.jnew("java/lang/String",.rnetlogo$savedworkingdir[2])),
			error = function(e) {
				e$printStackTrace()
      }
	)	
	
  if (gui) {
    .rnetlogo$startedGUI <- TRUE
    .rnetlogo$guiobj <- nl.obj
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

  # store NetLogo instance reference internally under the defined/submitted name
  assign(nl.obj, nlo, envir=.rnetlogo)
  .rnetlogo$objects <- c(.rnetlogo$objects, nl.obj)    
}

