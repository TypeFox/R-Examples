######################################################################
#
# zzz.R
#
# Edited by Zack Almquist
# Written by Carter T. Butts <buttsc@uci.edu>; based on an original by
# Carter T. Butts <buttsc@uci.edu>, David Hunter <dhunter@stat.psu.edu>,
# and Mark S. Handcock <handcock@u.washington.edu>.
# Last Modified 2/26/09
# Licensed under the GNU General Public License version 3 or later
#
# Part of the R/census package
#
# .First.lib is run when the package is loaded with library(UScensus2000tract)
#
######################################################################

.First.lib <- function(lib, pkg){
    if(R.version$major=="1"){
     ehelp <- help(package="UScensus2000tract")$info[[2]][[2]]
     cat(paste("'",ehelp[4],"'\n",
               "Version ",ehelp[2],
               " created on ",ehelp[3],".\n", sep=""))
    }else{
     ehelp <- help(package="UScensus2000tract")$info[[1]]
     cat(paste(substring(ehelp[3],first=16),"\n",
               "Version ",substring(ehelp[4],first=16),
               " created on ",
                substring(ehelp[5],first=16),".\n", sep=""))
    }
    cat(paste("copyright (c) 2009, Zack W. Almquist, University of California-Irvine\n",sep=""))
    cat('For citation information, type citation("UScensus2000tract").\n')
    cat('Type help(package="UScensus2000tract") to get started.\n')
}
