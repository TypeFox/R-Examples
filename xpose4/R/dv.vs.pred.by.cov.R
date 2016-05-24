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

"dv.vs.pred.by.cov" <-
  function(object,
           
           #xlb  = NULL,
           #ylb  = NULL,
           #onlyfirst=FALSE,
           #inclZeroWRES=FALSE,
           #subset=xsubset(object),
           #mirror=FALSE,
           #seed  = NULL,
           abline = c(0,1),
           #abllwd = object@Prefs@Graph.prefs$abllwd,
           #abllty = object@Prefs@Graph.prefs$abllty,
           #ablcol = object@Prefs@Graph.prefs$ablcol,
           #prompt = FALSE,
           smooth=TRUE,
           #main=NULL,
           #samp  = NULL,
           #max.plots.per.page=1,
           #scales=list(),
           #aspect = object@Prefs@Graph.prefs$aspect,#"fill"
           #mirror.aspect="fill",
           #pass.plot.list=FALSE,
           #x.cex=NULL,
           #y.cex=NULL,
           #main.cex=NULL,
           main="Default",
           ...) {

    
    ## Make sure we have the necessary variables defined in the 
    ## object.                                                  
    if(is.null(check.vars(c("dv","pred","covariates"),object))) {
      return(NULL)
    }

    
      ## create enpty list for plots
      number.of.plots <- 0
      for (i in xvardef("covariates", object)) {
          number.of.plots <- number.of.plots + 1
      }    
      plotList <- vector("list",number.of.plots)
      plot.num <- 0 # initialize plot number      
        
      ## loop
      for (i in xvardef("covariates", object)) {      


    
        xplot <- xpose.plot.default(xvardef("pred",object),
                                    xvardef("dv",object),
                                    #xlb = xlb,
                                    #ylb = ylb,
                                    abline=abline,
                                    #abllwd=abllwd,
                                    #scales=scales,
                                    #aspect=aspect,
                                    object,
                                    main=NULL,
                                    by=i,
                                    #subset=subset,
                                    smooth=smooth,
                                    #samp=samp,
                                    pass.plot.list=TRUE,
                                    ...)
                                  

        plot.num <- plot.num+1
        plotList[[plot.num]] <- xplot
      }

    default.plot.title <- paste(xlabel(xvardef("dv",object),object)," vs ",
                                xlabel(xvardef("pred",object),object),
                                sep="")
    plotTitle <- xpose.multiple.plot.title(object=object,
                                           plot.text = default.plot.title,
                                           main=main,
                                           ...)
    obj <- xpose.multiple.plot(plotList,plotTitle,...)
    return(obj)

  }
