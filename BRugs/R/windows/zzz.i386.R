".onLoad.i386" <- function(lib, pkg){
    ob <- findOpenBUGS()
    loadOpenBUGS(ob$dir)
    msg <- paste("Welcome to BRugs connected to OpenBUGS")
    if (!is.na(ob$version))
        msg <- paste(msg, "version", ob$version)
    else msg <- paste(msg, "in directory", ob$dir)
    packageStartupMessage(msg)
}

".onUnload.i386" <- function(libpath){
    if(is.loaded("CmdInterpreter")) {
        libname <- paste(options()$OpenBUGS, "libOpenBUGS.dll", sep="/")
        dyn.unload(libname)
    }
}

## Load OpenBUGS from specified location
loadOpenBUGS <- function(dir)  {
    libname <- paste(dir, "libOpenBUGS.dll", sep="/")
    if (!file.exists(libname))  {
        warning("Shared library \"libOpenBUGS.dll\" not found in ", dir)
        return(FALSE)
    }
    options(OpenBUGS = dir)
    dyn.load(libname)
    len <- nchar(dir)
    .C("SetWorkingDir", as.character(dir), len, PACKAGE="libOpenBUGS")
    ## Set temporary dir for "buffer.txt" output
    tempDir <- gsub("\\\\", "/", tempdir())
    .C("SetTempDir", as.character(tempDir), nchar(tempDir), PACKAGE="libOpenBUGS")
    command <- "BugsMappers.SetDest(2)"
    .CmdInterpreter(command)
    if(is.null(getOption("BRugsVerbose")))
        options("BRugsVerbose" = TRUE)
    .initGlobals()
    options(OpenBUGSExamples = paste(dir, "Examples", sep="/"))
    invisible()
}
