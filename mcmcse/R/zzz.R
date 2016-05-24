######################################################################
#
# zzz.R
#
# Written by John Hughes <hughesj@umn.edu>.
#
# Last Modified 08/04/12
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/mcmcse package
#
# .onLoad is run when the package is loaded with library(mcmcse)
#
######################################################################

#' @import utils

.onAttach = function(libname, pkgname)
{
    # library.dynam("batchmeans", package = pkgname, lib.loc = libname)
    temp = packageDescription("mcmcse")
    msg = paste(temp$Package, ": ", temp$Title, "\n", "Version ", temp$Version,
                " created on ", temp$Date, ".\n", sep = "")
    msg = paste(msg,
"copyright (c) 2012, James M. Flegal, University of California,Riverside\n",
"                    John Hughes, University of Minnesota\n",
"                    Dootika Vats, University of Minnesota\n", sep = "")
    msg = paste(msg, 'For citation information, type citation("mcmcse").\n')
    msg = paste(msg, 'Type help("mcmcse-package") to get started.\n')
    packageStartupMessage(msg)
}
