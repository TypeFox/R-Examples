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

"ind.plots.wres.hist" <-
  function(object,
           main = "Default",
           wres="wres",
           #xlb  = NULL,
                                        # ylb  = xlabel(xvardef("dv",object),object),
           ylb = NULL,
           layout=c(4,4),
           inclZeroWRES=FALSE,
           subset=xsubset(object),
           scales=list(cex=0.7,tck=0.5),
           aspect="fill",
           force.by.factor=TRUE,
           ids=F,
           as.table=TRUE,
           hicol = object@Prefs@Graph.prefs$hicol,
           hilty = object@Prefs@Graph.prefs$hilty,
           hilwd = object@Prefs@Graph.prefs$hilwd,
           hidcol = object@Prefs@Graph.prefs$hidcol,
           hidlty = object@Prefs@Graph.prefs$hidlty,
           hidlwd = object@Prefs@Graph.prefs$hidlwd,
           hiborder = object@Prefs@Graph.prefs$hiborder,
           prompt = FALSE,
           mirror=NULL,
           main.cex=0.9,
           max.plots.per.page=1,
                                        #lty=c(0,1,1),
                                        #pch=c(21,32,32),
                                        #type="o",
                                        #col=c(1,object@Prefs@Graph.prefs$smcol,object@Prefs@Graph.prefs$lmcol),
                                        #lwd=1,
           ...) {

    ## Make sure we have the necessary variables defined in the ##
    ## object.                                                  ##
    if(is.null(check.vars(c("id",wres),object))) {
      return(NULL)
    }
    

    ## check for mirror
    if(!is.null(mirror)){
      cat("Mirror not currently implemented for individual plots\n")
      return()
    }

    data <- Data(object,inclZeroWRES,subset=subset)
    
    ## Bin them
    list.id   <- unique(data[[xvardef("id",object)]])
    length.id <- length(list.id)
    plots.per.page <- layout[1] * layout[2]
    plots.cur <- 0
    pages <- 1
    page.breaks <- c(0)
    old.obj <- object
    new.obj <- object
    new.obj@Data <- NULL
    
    for (i in list.id) {
      plots.cur <- plots.cur + 1
      if (plots.cur == plots.per.page) {
        pages <- pages + 1
        plots.cur <- 0
        page.breaks <- c(page.breaks, i)
      }
    }
    if (max(page.breaks) < max(list.id)) {
      page.breaks <- c(page.breaks, max(list.id))
    }
    data$bin <- cut(data$ID, page.breaks, include.lowest=T)
    id.levels <- levels(data$bin)
    

    plot.num <- 0
    plotList <- vector("list",length(id.levels))     
    for (i in id.levels) {    ## start loop
      
      new.obj@Data <- data[data$bin==i,] #subset(data, bin == i)

      ## Set up the data ##
      ## nobj <- new("xpose.data",
      ##             Runno=object@Runno,
      ##             Data = NULL 
      ##             )
      ## Data(nobj) <- Data(new.obj,inclZeroWRES=inclZeroWRES,
      ##                    subset=subset)

      if(is.null(xvardef(wres,object))){
        plotvar <- wres
      }else{
        plotvar <- xvardef(wres,object)
      }

      ## Fix any main and/or axis titles
      default.plot.title <- paste("Individual plots of", plotvar, sep=" ")
      plotTitle <- xpose.multiple.plot.title(object=object,
                                             plot.text = default.plot.title,
                                             main=main,
                                             ...)
      
                                        # Set y axis title
##       if (is.null(xlb)) {
##         xlb <- xlabel(xvardef("wres",object),object)
##       }
##       if (is.null(ylb)) {
##         ylb <- "Proportion"
##       }
      
      
      xplot <- xpose.plot.histogram(plotvar,#xvardef("wres",nobj),
                                    new.obj,
                                    #xlb = xlb,
                                    #ylb = ylb,
                                    by=xvardef("id",new.obj),
                                    main=plotTitle,
                                        #group="ind",
                                    layout=layout,
                                    scales=scales,
                                    aspect=aspect,
                                    xvar = plotvar,#xvardef("wres",object),
                                    force.by.factor=force.by.factor,
                                    ids=ids,
                                    subset=subset,
                                    as.table=as.table,
                                    hicol = hicol,
                                    hilty = hilty,
                                    hilwd = hilwd,
                                    hidcol = hidcol,
                                    hidlty = hidlty,
                                    hidlwd = hidlwd,
                                    hiborder = hiborder,
                                    main.cex=main.cex,
                                    ...)
      plot.num <- plot.num+1
      plotList[[plot.num]] <- xplot
      

    }

    obj <- xpose.multiple.plot(plotList,max.plots.per.page=max.plots.per.page,plotTitle=NULL,prompt=prompt,...)
    return(obj)


  }
