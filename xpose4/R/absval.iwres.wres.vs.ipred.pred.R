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

"absval.iwres.wres.vs.ipred.pred" <-
  function(object,
           
           ##aspect="fill",
           main="Default",
           ...) {

  
    if(is.null(check.vars(c("pred","wres","iwres","ipred"),
                          object,silent=FALSE))) {
      return()
    }

    num.of.plots <- 2
    plotList <- vector("list",num.of.plots)

    plot1 <- absval.wres.vs.pred(object,main=NULL,
                              ##aspect=aspect,
                              pass.plot.list=TRUE,
                              ...)
    plot2 <- absval.iwres.vs.ipred(object,main=NULL,
                                ##aspect=aspect,
                                pass.plot.list=TRUE,
                                ...)

    plotList[[1]] <- plot1
    plotList[[2]] <- plot2

    default.plot.title <- "Weighted residuals vs. Predictions"
    
    plotTitle <- xpose.multiple.plot.title(object=object,
                                           plot.text = default.plot.title,
                                           main=main,
                                           ...)
    obj <- xpose.multiple.plot(plotList,plotTitle,...)
    return(obj)

}
