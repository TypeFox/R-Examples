# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

"absval.cwres.vs.cov.bw" <-
  function(object,
           
           xlb = "|CWRES|",
           main="Default",
           ...) {

    if (is.null(xvardef("covariates", object))) {
        return(cat("Covariates are not properly set in the database!\n"))
    }
    
    if(is.null(xvardef("cwres",object))) {
        return(cat("There are no CWRES defined in the database!\n"))    
    }


    ## create list for plots
    number.of.plots <- 0
    for (i in xvardef("covariates", object)) {
      number.of.plots <- number.of.plots + 1
    }
    plotList <- vector("list",number.of.plots)
    
    
    ## loop through covariates
    plot.num <- 0 # initialize plot number
    for (i in xvardef("covariates", object)) {
    
      xplot <- xpose.plot.bw(xvardef("cwres",object),
                             i,
                             object,
                             main=NULL,
                             xlb=xlb,
                             binvar = i,
                             funx="abs",
                             pass.plot.list = TRUE,
                             ...)
      

      plot.num <- plot.num+1
      plotList[[plot.num]] <- xplot
    }
    
    default.plot.title <- paste("|",xlabel(xvardef("cwres",object),object),"| vs ",
                                "Covariates", sep="")
    
    plotTitle <- xpose.multiple.plot.title(object=object,
                                           plot.text = default.plot.title,
                                           main=main,
                                           ...)
    obj <- xpose.multiple.plot(plotList,plotTitle,...)
    return(obj)

  }
