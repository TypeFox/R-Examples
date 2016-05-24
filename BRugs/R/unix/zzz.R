if (is.R()){

    ".onLoad" <- function(lib, pkg){
        ## TODO any need for these to be user specifiable?
        options("BRugsTmpdir" = gsub("\\\\", "/", tempdir()))
        options("BRugsExtFile" = paste(basename(tempfile()), ".bug", sep=""))
        options(OpenBUGS = "/usr/local/lib/OpenBUGS/lib")
        options(OpenBUGSdoc = "/usr/local/lib/OpenBUGS/doc")
        options(OpenBUGSExamples = paste(options()$OpenBUGSdoc, "Examples", sep="/"))

        if(is.null(getOption("BRugsVerbose")))
            options("BRugsVerbose" = TRUE)
        .initGlobals()
        ver <- system("echo \"modelQuit()\" | /usr/local/lib/OpenBUGS/lib/../bin/OpenBUGS", intern=TRUE)
        ver <- sub("OpenBUGS version (([0-9]\\.)+[0-9]).+","\\1",ver[1])
        packageStartupMessage("Welcome to BRugs connected to OpenBUGS version ", ver)
    }

    ".onUnload" <- function(libpath){
    }

    ## Windows-only
    loadOpenBUGS <- function(dir)  {
    }
}
