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

"ranpar.qq" <-
  function(object,
           onlyfirst=TRUE,
           main="Default",
           
           ...) {
    
    ## is everything in place?
    if(is.null(xvardef("ranpar",object))) {
      return(cat("There are no ETAs in the database!\n"))    
    }
    
    ## create enpty list for plots
    number.of.plots <- 0
    for (i in xvardef("ranpar", object)) {
      if(!is.factor(object@Data[[i]])){
        number.of.plots <- number.of.plots + 1
      }
    }    
    plotList <- vector("list",number.of.plots)
    plot.num <- 0 # initialize plot number

    ## loop (ranpar)
    for (i in xvardef("ranpar", object)) {      
      
      if(!is.factor(object@Data[[i]])){

        xplot <- xpose.plot.qq(i,
                               object,
                               main=NULL,
                               onlyfirst=onlyfirst,
                               pass.plot.list=TRUE,
                               ...)
        
        
        plot.num <- plot.num+1
        plotList[[plot.num]] <- xplot
      }
    }

    
    default.plot.title <- "Distribution of random parameters"
    plotTitle <- xpose.multiple.plot.title(object=object,
                                           plot.text = default.plot.title,
                                           main=main,
                                           ...)
    obj <- xpose.multiple.plot(plotList,plotTitle,...)
    return(obj)    
    
  }
