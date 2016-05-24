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

"ind.plots" <-
  function(object,
           y.vals = c(
             xvardef("dv",new.obj),
             xvardef("ipred",new.obj),
             xvardef("pred",new.obj)
           ),
           x.vals = xvardef("idv",new.obj),
           id.vals = xvardef("id",new.obj),
           key.text = y.vals,
           main = "Default",
           key="Default",
           xlb = xlabel(xvardef("idv",object),object),
           ylb = NULL,
           layout = c(4,4),
           inclZeroWRES = FALSE,
           subset = xsubset(object),
           type = "o",
           grid=FALSE,
           col = c(1,2,4),
           lty = c(0,1,3),
           lwd = c(1,1,1),
           pch = c(21,32,32),
           cex = c(0.7,0.7,0.7),
           fill=c("lightgrey",0,0),
           prompt = FALSE,
           mirror=NULL,
           main.cex=0.9,
           max.plots.per.page=1,
           pch.ip.sp=c(21,19,18),
           cex.ip.sp=c(0.7,0.4,0.4),
           y.vals.subset=NULL,
           ...){
    
    ## check for mirror
    if(!is.null(mirror)){
      cat("Mirror not currently implemented for individual plots\n")
      return()
    }
    
    ## Make sure we have the necessary variables defined in the ##
    ## object.                                                  ##
    if(is.null(check.vars(c("id","idv","dv","ipred","pred"),object))) {
      return(NULL)
    }
    
    
    ## Bin them
    old.obj <- object 
    old.obj@Data <- Data(object,inclZeroWRES,onlyfirst=FALSE,subset=subset)
    
    if(any(is.null(old.obj@Data))){
      return("The subset expression produces no data.")
    }
    
    list.id   <- sort(unique(old.obj@Data[[xvardef("id",object)]]))
    length.id <- length(list.id)
    plots.per.page <- layout[1] * layout[2]
    plots.cur <- 0
    pages <- 1
    page.breaks <- c(0)
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
    id.levels <- levels(cut(old.obj@Data[[xvardef("id",object)]], page.breaks, include.lowest=T))
    old.obj@Data$bin <- cut(old.obj@Data[[xvardef("id",object)]], page.breaks, include.lowest=T)
    
    plot.num <- 0
    plotList <- vector("list",length(id.levels))
    for (i in id.levels) { 
      
      new.obj@Data <- old.obj@Data[old.obj@Data$bin==i,] #==subset(old.obj@Data, bin == i)
      
      ## Set up the data ##      
      ## Figure out what variables we have defined
      select <- y.vals
      key.text <- key.text
      
      numpans <- length(select)
      
      nobj <- new("xpose.data",
                  Runno=object@Runno,
                  Data = NULL)
      
      Data(nobj,keep.structure=T) <- xpose.stack(Data(new.obj,inclZeroWRES=inclZeroWRES,
                                                      subset=subset),
                                                 new.obj,
                                                 select=select,
                                                 rep=c(
                                                   x.vals,
                                                   id.vals
                                                 ),
                                                 subset=y.vals.subset)
      
      ## Fix any main and/or axis titles
      default.plot.title <- "Individual plots"
      plotTitle <- xpose.multiple.plot.title(object=object,
                                             plot.text = default.plot.title,
                                             main=main,
                                             subset=subset,
                                             ...)
      
      ## Set y axis title
      if(is.null(ylb)) {
        ylb <- "Observations / Predictions"
      }
      
      
      if(!all(is.na(match(key,"Default")))) {
        key=list(text=list(key.text),columns=numpans,
                 #points=list(pch=pch,cex=cex),
                 #rectangles=list(pch=pch,cex=cex),
                 lines=list(lty=lty,pch=pch,cex=cex,col=col,lwd=lwd,size=3,type=type,fill=fill),
                 between=1)
      }
      
      xplot <- xpose.plot.default(x.vals,
                                  "values",
                                  nobj,
                                  by=id.vals,
                                  groups="ind",
                                  iplot=TRUE,
                                  layout=layout,
                                  main=plotTitle,
                                  inclZeroWRES=TRUE,
                                  xlb = xlb,
                                  ylb = ylb,
                                  lty=lty,
                                  col=col,
                                  cex=cex,
                                  grid=grid,
                                  pch=pch,
                                  type=type,
                                  lwd=lwd,
                                  fill=fill,
                                  force.by.factor=TRUE,
                                  ids=F,
                                  as.table=TRUE,
                                  key = key,
                                  subset=subset,
                                  main.cex=main.cex,
                                  force.x.continuous=TRUE,
                                  pch.ip.sp=pch.ip.sp,
                                  cex.ip.sp=cex.ip.sp,
                                  ...)
      plot.num <- plot.num+1
      plotList[[plot.num]] <- xplot
      
      
    }
    obj <- xpose.multiple.plot(plotList,plotTitle=NULL,max.plots.per.page=max.plots.per.page,prompt=prompt,...)
    return(obj)
    
  }
