######################################################################
#
# zzz.R
#
# Written by John Hughes <hughesj@umn.edu>.
#
# Last Modified 01/16/15
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/ngspatial package
#
# .onLoad is run when the package is loaded with library(ngspatial)
#
######################################################################

#' @import batchmeans Rcpp utils

.onLoad = function(libname, pkgname)
{
    loadRcppModules(direct = FALSE)
}

.onAttach = function(libname, pkgname)
{
    temp = packageDescription("ngspatial")
    msg = paste(temp$Package, ": ", temp$Title, "\n", "Version ", temp$Version,
                " created on ", temp$Date, ".\n", sep = "")
    msg = paste(msg, "copyright (c) 2013-15, John Hughes, University of Minnesota\n",
"                       Xiaohui Cui, University of Minnesota\n",
                sep = "")
    msg = paste(msg, 'For citation information, type citation("ngspatial").\n', sep = "")
    msg = paste(msg, 'Type help(package = ngspatial) to get started.\n', sep = "")
    packageStartupMessage(msg)
}
