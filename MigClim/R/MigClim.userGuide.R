#
# Displays the PDF user guide to the MigClim R package.
#
MigClim.userGuide <- function()
{
    UGpath <- find.package("MigClim")
    if(.Platform$OS.type=="windows") eval(parse(text=paste("system('open \"", UGpath, "/doc/MigClim_userGuide.pdf\"')", sep="")))
    if(.Platform$OS.type=="unix"){
        cat("If the user guide does not open automatically, it can be found in the following directory: ", UGpath, "/doc/MigClim_userGuide.pdf \n", sep="")
        
        if(Sys.info()["sysname"] =="Linux"){
            system(paste("evince ", UGpath, "/doc/MigClim_userGuide.pdf", sep=""))  # works with linux system that have evince installed.
        } else system(paste("open ", UGpath, "/doc/MigClim_userGuide.pdf", sep="")) # should work for Mac OS but not sure.  
    }
}


