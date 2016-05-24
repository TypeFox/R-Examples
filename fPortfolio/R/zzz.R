
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################


.onAttach <- 
function(libname, pkgname)
{
  # do whatever needs to be done when the package is loaded
  # some people use it to bombard users with 
  # messages using 
  
  packageStartupMessage( "\n" )
  packageStartupMessage( "Rmetrics Package fPortfolio" ) 
  packageStartupMessage( "Portfolio Optimization" )
  packageStartupMessage( "Copyright (C) 2005-2014 Rmetrics Association Zurich" )  
  packageStartupMessage( "Educational Software for Financial Engineering and Computational Science" ) 
  packageStartupMessage( "Rmetrics is free software and comes with ABSOLUTELY NO WARRANTY." ) 
  packageStartupMessage( "https://www.rmetrics.org --- Mail to: info@rmetrics.org" ) 
}


###############################################################################


.onLoad <- 
  function(libname, pkgname)
{
    if(!is.numeric(timeDate::getRmetricsOptions("length.print"))) 
        timeDate::setRmetricsOptions(length.print = 5)
    
    timeDate::setRmetricsOptions(.x.save = NA)
    
    eval(attach <- function(what) 
      base::attach(what, warn.conflicts=FALSE), envir=.GlobalEnv)

}
    # Startup Mesage and Desription:
    # MSG <- if(getRversion() >= "2.5") packageStartupMessage else message
    # dsc <- packageDescription(pkg)
    # if(interactive() || getOption("verbose")) { 
    #    title <- paste(strsplit(dsc$Title, split = "-")[1:2])
    #    MSG(paste(
    #        "\nPackage ", pkg, " (", dsc$Version, ") loaded.\n",
    #        dsc$Title, "\n", 
    #        dsc$Copyright, ", ", dsc$License, "\n", 
    #        dsc$Author, "\n", 
    #        dsc$URL, "\n", sep="")) 
    # }


###############################################################################


