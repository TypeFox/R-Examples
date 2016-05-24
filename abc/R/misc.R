######################################################################
#
# misc.R
#
# copyright (c) 2011-05-30, Katalin Csillery, Olivier Francois and
# Michael GB Blum with some initial code from Mark Beaumont
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
# 
# Part of the R/abc package
# Contains: normalise,namesWarningFilter
#
######################################################################


normalise <- function(x,y){
  if(mad(y) == 0)
    return (x)
  else
    return (x/mad(y))
}


namesWarningFilter <- function(x){
    if( any( grepl( "No parameter names are given", x) ) ) invokeRestart( "muffleWarning" )
    if( any( grepl( "No summary statistics names are given", x) ) ) invokeRestart( "muffleWarning" )
}
