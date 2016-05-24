if (is.R()){

    ".onLoad.x64" <- function(lib, pkg){
        ## TODO any need for these to be user specifiable?
        options("BRugsTmpdir" = gsub("\\\\", "/", tempdir()))
        options("BRugsExtFile" = paste(basename(tempfile()), ".bug", sep=""))
        ob <- findOpenBUGS()
        options(OpenBUGS = ob$dir)
        options(OpenBUGSdoc = ob$dir)
        options(OpenBUGSExamples = paste(ob$dir, "Examples", sep="/"))

        if(is.null(getOption("BRugsVerbose")))
            options("BRugsVerbose" = TRUE)
        .initGlobals()
        msg <- paste("Welcome to BRugs connected to OpenBUGS")
        if (!is.na(ob$version))
            msg <- paste(msg, "version", ob$version)
        else msg <- paste(msg, "in directory", ob$dir)
        packageStartupMessage(msg)
        pathtoBUGS <- gsub("/", "\\", ob$dir)
        oldpath <- Sys.getenv("PATH")
        if(!length(grep(pathtoBUGS, oldpath, fixed=TRUE)))
            Sys.setenv(PATH=paste(oldpath, pathtoBUGS, sep=";"))
    }

    ".onUnload.x64" <- function(libpath){
    }

}
