## generic.R
## Author          : Claus Dethlefsen
## Created On      : Mon Nov 19 20:48:24 2001
## Last Modified By: Claus Dethlefsen
## Last Modified On: Tue Nov 01 13:46:57 2011
## Update Count    : 109
## Status          : Unknown, Use with caution!
###############################################################################
##
##    Copyright (C) 2002  Susanne Gammelgaard Bottcher, Claus Dethlefsen
##
##    This program is free software; you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program; if not, write to the Free Software
##    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
######################################################################



printline <- function(s="-",n=60) cat(rep(s,n),"\n",sep="")

#.First.lib <- function(lib, pkg)
#{
#    require(methods)
#    require(dynamicGraph)
#    library.dynam("deal", package = pkg, lib.loc = lib)

#      if((R.version$major == 1) && (as.numeric(R.version$minor) < 9))
#        packageDescription <- package.description

#    cat("\n")
#    cat("-------------------------------------------------------------\n")
#    cat(packageDescription("deal", lib = lib, field="Title"))
#    cat("\n")
#    ver <- packageDescription("deal", lib = lib, field="Version")
#    maint<- packageDescription("deal", lib = lib, field="Maintainer")
#    built<- packageDescription("deal", lib = lib, field="Built")
#    URL  <- packageDescription("deal", lib = lib, field="URL")
#    cat(paste("deal, version", ver,  "is now loaded\n"))
#    cat("Copyright (C) 2002-2007, Susanne G. Bottcher and Claus Dethlefsen\n")
#    cat("Maintained by",maint,"\n")
#    cat("Webpage:",URL,"\n")
#    cat("\nBuilt:",built,"\n")
#    cat("-------------------------------------------------------------\n")
#    cat("\n")

#    require(methods)
#    .load.deal.networkclass()
#    .load.dynamicgraph()
#  return(invisible(0))
#}

.onAttach <- function (lib, pkg) 
{
#    require(methods)
#    .load.deal.networkclass()
# library.dynam("deal", package = pkg, lib.loc = lib)
  }

.onLoad <- function (lib, pkg) 
{
#    require(methods)
#    .load.deal.networkclass()
library.dynam("deal", package = pkg, lib.loc = lib)
}


#.Last.lib <- function(lib) {
#  cat("Thank you for using deal\n")
#  return(invisible(0))
#}

