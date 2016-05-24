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

"addit.gof" <-
  function(object,
           type="p",
           title.size=0.02,
           title.just=c("center","top"),
           
           main="Default",
           force.wres=FALSE,
           ...) {

    if(is.null(check.vars(c("dv","pred","ipred","iwres","idv"),
                          object,silent=FALSE))) {
      return()
    }

    use.cwres=TRUE
    if(force.wres){
      use.cwres=FALSE
      if(is.null(check.vars(c("wres"),object,silent=FALSE))) return()
    } else {
      if(is.null(check.vars(c("cwres"),object,silent=TRUE))) {
        use.cwres=FALSE
        if(is.null(check.vars(c("wres"),object,silent=FALSE))) return()
      }
    }
    
    ## create empty list for plots
    num.of.plots <- 4
    plotList <- vector("list",length(num.of.plots))

    if(use.cwres){
      xplot1 <- cwres.vs.pred(object,
                              main=NULL,
                              type=type,
                              pass.plot.list=TRUE,
                              ...)
      
      xplot3 <- cwres.vs.idv(object,
                            main=NULL,
                            type="n",
                            ids=F,
                            pass.plot.list=TRUE,
                            ...)

      xplot4 <- cwres.vs.idv(object,
                            main=NULL,
                            logx=T,
                            type=type,
                            pass.plot.list=TRUE,
                            ...)
      

    } else {
      xplot1 <- wres.vs.pred(object,
                             main=NULL,
                             type=type,
                             pass.plot.list=TRUE,
                             ...)

      xplot3 <- wres.vs.idv(object,
                            main=NULL,
                            type="n",
                            ids=F,
                            pass.plot.list=TRUE,
                            ...)

      xplot4 <- wres.vs.idv(object,
                            main=NULL,
                            logx=T,
                            type=type,
                            pass.plot.list=TRUE,
                            ...)
      
    }
    
    xplot2 <- absval.iwres.vs.idv(object,
                                  main=NULL,
                                  ylb=paste("Absolute value of iWRES"),
                                  pass.plot.list=TRUE,
                                  ...)

    
    plotList[[1]] <- xplot1
    plotList[[2]] <- xplot2
    plotList[[3]] <- xplot3
    plotList[[4]] <- xplot4      


    default.plot.title <- "Additional goodness of fit plots"
    plotTitle <- xpose.multiple.plot.title(object=object,
                                           plot.text = default.plot.title,main=main,
                                           ...)

    
    obj <- xpose.multiple.plot(plotList,plotTitle,...)
    return(obj)

  }

