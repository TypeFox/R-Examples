.First.lib <- function(lib, pkg){
    if(R.version$major=="1"){
     ehelp <- help(package="yacca")$info[[2]][[2]]
     cat(paste("'",ehelp[4],"'\n",
               "Version ",ehelp[2],
               " created on ",ehelp[3],".\n", sep=""))
    }else{
     ehelp <- help(package="yacca")$info[[1]]
     cat(paste(substring(ehelp[3],first=16),"\n",
               "Version ",substring(ehelp[4],first=16),
               " created on ",
                substring(ehelp[5],first=16),".\n", sep=""))
    }
    cat("Copyright (c) 2008, Carter T. Butts, University of California-Irvine\n")
    cat('For citation information, type citation("yacca").\n')
    cat('Type help("yacca-package") to get started.\n')
}
