######################################################################
#
# zzz.r
#
# copyright (c) 2004, Mark S. Handcock, University of Washington
# written June 2004
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/glmc package
#
# .First.lib is run when the package is loaded with library(glmc)
#
######################################################################

#.First.lib <- function(lib, pkg){ 
#   .onAttach<-function(lib, pkg){
#   library.dynam("glmc", pkg, lib)
#    ehelp <- library(help="glmc",lib.loc=NULL,character.only=TRUE)$info[[1]]
#    cat(paste(substring(ehelp[4],first=16),"\n",
#              "Version ",substring(ehelp[2],first=16),
#              " created on ", substring(ehelp[3],first=16),".\n", sep=""))
#    cat(paste("'",ehelp[4],"'\n",
#              "Version ",ehelp[2],
#              " created on ",ehelp[3],".\n", sep=""))
#    cat(paste("copyright (c) 2004, Mark S. Handcock, University of Washington\n",sep=""))
#    cat('To cite, see citation("glmc")\n')
#    cat('Type help(package="glmc") to get started.\n')
#    require(emplik)
#}

.onAttach<-function(lib, pkg){
ehelp <- library(help="glmc",lib.loc=NULL,character.only=TRUE)$info[[1]]
packageStartupMessage(paste(substring(ehelp[4],first=16),"\n",
              "Version ",substring(ehelp[2],first=16),
              " created on ", substring(ehelp[3],first=16),".\n", sep=""))
    packageStartupMessage(paste("'",ehelp[4],"'\n",
              "Version ",ehelp[2],
              " created on ",ehelp[3],".\n", sep=""))
    packageStartupMessage(paste("copyright (c) 2004, Mark S. Handcock, University of Washington\n",sep=""))
    packageStartupMessage('To cite, see citation("glmc")\n')
    packageStartupMessage('Type help(package="glmc") to get started.\n')
}

