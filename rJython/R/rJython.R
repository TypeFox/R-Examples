#########################################################
# CGB, 20100627
#########################################################

rJython <- function( jython.jar = NULL, modules = NULL ){

    stopifnot(require(rJava))

    # Aux function
    
    system.file. <- function(...) {
        s <- system.file(...)
        if (.Platform$OS == "windows") gsub("\\", "/", s, fixed = TRUE) else s
    }

    # Looking for jython jar 

    if( is.null( jython.jar ) ) 
        jython.jar <- Sys.getenv("RJYTHON_JYTHON")
    if (is.null(jython.jar) || jython.jar == "")
        jython.jar <- system.file.("jython.jar", package = "rJython")

    # Setting the temp dir
   
    tmp.dir.parm <- paste( "-Djava.io.tmpdir=", tempdir() )

    # Starting JVM

    .jinit(jython.jar, parameters = tmp.dir.parm )
    rJython <- .jnew("org.python.util.PythonInterpreter")

    # Adding required python modules to the interpreter

    rJython$exec( "import sys" )

    if (is.character(modules)) modules <- as.list(modules)
    modules <- c( modules, list( system.file.( package = "rJython" ) ) )
    modules <- lapply( modules, function( module ) paste( "sys.path.append(", module, ");", sep = '"' ) )
    lapply( modules, rJython$exec )
    
    rJython$exec( "import simplejson as json" )

    rJython
}

