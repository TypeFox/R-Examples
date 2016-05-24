## See unix/zzz.R, windows/zzz.R for platform specific .onLoad functions

if (is.R()){

    .initGlobals <- function(){
        options("BRugsSamplesBeg" = 1)
        options("BRugsSamplesEnd" = 10000000)
        options("BRugsSamplesFirstChain" = 1)
        options("BRugsSamplesLastChain" = 1)
        options("BRugsSamplesThin" = 1)
        options("BRugsSamplesVariable" = "*")
        options("BRugsNextChain" = 1) # index of chain which needs to be initialized next
        options("BRugsPrec" = 4)
    }

    ## Overwriting new (from R-2.6.0) sQuote (for typing human readable text) in R within the BRugs Namespace!
    ## we cannot use sQuote that uses fancy quotes!
    sQuote <- function(x) paste("'", x, "'", sep="")


} else {  # ends if (is.R())

    ".First.lib" <- function(lib.loc, section)
    {
        dyn.open(system.file("OpenBUGS", "brugs.dll", package="BRugs"))
        ## sets path / file variables and initializes subsystems
        root <- file.path(system.file("OpenBUGS", package="BRugs"))
        len <- nchar(root)
        tempDir <- gsub("\\\\", "/", tempdir())
        .C("SetRoot", as.character(root), len)
        .C("SetTempDir", as.character(tempDir), nchar(tempDir))
        command <- "BugsMappers.SetDest(2)"
        .C("CmdInterpreter", as.character(command), nchar(command), integer(1))
        if(is.null(getOption("BRugsVerbose")))
            options("BRugsVerbose" = TRUE)
        invisible()
    }

    .tempDir <- getwd()

    tempdir <- function(){ .tempDir }

}  # ends else
