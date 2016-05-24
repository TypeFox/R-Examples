.onAttach <-
function(libname, pkgname)
{
    mylib <- dirname(system.file(package = "tseries"))
    ver <- packageDescription("tseries", lib.loc = mylib)["Version"]
    txt <- c("\n",
             paste(sQuote("tseries"), "version:", ver),
             "\n",
             paste(sQuote("tseries"),
                   "is a package for time series analysis",
                   "and computational finance."),
             "\n",
             paste("See",
                   sQuote("library(help=\"tseries\")"),
                   "for details."),
             "\n")
    if(interactive() || getOption("verbose"))
        packageStartupMessage(paste(strwrap(txt,
                                            indent = 4,
                                            exdent = 4),
                                    collapse = "\n"))
                                    
}
