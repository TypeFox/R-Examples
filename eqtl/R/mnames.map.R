#####################################################################
#
# mnames.map.R
#
# copyright (c) 2008-3, Ahmid A Khalili
# 
# last modified Jul, 2008
# first written Mar, 2008
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/eqtl package
# Contains: mnames.map
#
######################################################################

######################################################################
#
# mnames.map: List all markers
#
######################################################################

`mnames.map` <-
function(cross)
{

require(qtl)

if ( class(cross)[2] != "cross" ) stop("peak should have class \"cross\"\n")

mnames <- NA
for ( i in seq(names(cross$geno)) ) mnames <- c(mnames,names(cross$geno[[i]]$map))

try(return(mnames[-1]),silent=FALSE)

}

