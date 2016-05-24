######################################################################
#
# zzz.R
#
# Written by John Hughes <hughesj@umn.edu>.
#
# Last Modified 01/16/15
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/batchmeans package
#
# .onLoad is run when the package is loaded with library(batchmeans)
#
######################################################################

#' @import utils

.onAttach = function(libname, pkgname)
{
    # library.dynam("batchmeans", package = pkgname, lib.loc = libname)
    temp = packageDescription("batchmeans")
    msg = paste(temp$Package, ": ", temp$Title, "\n", "Version ", temp$Version,
                " created on ", temp$Date, ".\n", sep = "")
    msg = paste(msg,
"copyright (c) 2012-15, Murali Haran, Penn State University\n",
"                    John Hughes, University of Minnesota\n", sep = "")
    msg = paste(msg, 'For citation information, type citation("batchmeans").\n', sep = "")
    msg = paste(msg, 'Type help(package = batchmeans) to get started.\n', sep = "")
    packageStartupMessage(msg)
}
