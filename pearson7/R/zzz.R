######################################################################
#
# zzz.R
#
# Written by John Hughes <hughesj@umn.edu>.
#
# Last Modified 01/16/15
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/pearson7 package
#
# .onLoad is run when the package is loaded with library(pearson7)
#
######################################################################

#' @import utils

.onAttach = function(libname, pkgname)
{
    temp = packageDescription("pearson7")
    msg = paste(temp$Package, ": ", temp$Title, "\n", "Version ", temp$Version,
                " created on ", temp$Date, ".\n", sep = "")
    msg = paste(msg, "copyright (c) 2013-15, John Hughes, University of Minnesota\n",
                sep = "")
    msg = paste(msg, 'For citation information, type citation("pearson7").\n', sep = "")
    msg = paste(msg, 'Type help(package = pearson7) to get started.\n', sep = "")
    packageStartupMessage(msg)
}

