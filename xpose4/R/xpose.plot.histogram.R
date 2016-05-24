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

xpose.plot.histogram <-
  function(x,object,
           inclZeroWRES = FALSE,
           onlyfirst    = FALSE,
           samp         = NULL,
           type         = "density",
           aspect       = object@Prefs@Graph.prefs$aspect,
           scales       = list(),
           
           ## Conditioning settings
           by           = object@Prefs@Graph.prefs$condvar,
           force.by.factor = FALSE,
           ordby     = object@Prefs@Graph.prefs$ordby,
           byordfun  = object@Prefs@Graph.prefs$byordfun,
           shingnum  = object@Prefs@Graph.prefs$shingnum,
           shingol   = object@Prefs@Graph.prefs$shingol,
           strip = function(...)
           strip.default(...,strip.names=c(TRUE,TRUE)),
           ##par.strip.text=trellis.par.get("add.text"),
           
           ## Subset stuff
           subset       = xsubset(object),
           ## Axes and titles
           ##main         = NULL,
           main  = xpose.create.title.hist(x,object,subset,...),
           xlb          = NULL,
           ylb          = "Density", # this should be dependent on type
           
           ## Colors and stuff
           hicol = object@Prefs@Graph.prefs$hicol,
           hilty = object@Prefs@Graph.prefs$hilty,
           hilwd = object@Prefs@Graph.prefs$hilwd,
           hidcol = object@Prefs@Graph.prefs$hidcol,
           hidlty = object@Prefs@Graph.prefs$hidlty,
           hidlwd = object@Prefs@Graph.prefs$hidlwd,
           hiborder = object@Prefs@Graph.prefs$hiborder,

           ## mirror stuff
           mirror       = FALSE,
           max.plots.per.page=4,
           mirror.aspect="fill",
           pass.plot.list=FALSE,
           x.cex=NULL,
           y.cex=NULL,
           main.cex=NULL,
           mirror.internal=list(strip.missing=missing(strip)),
           ...) {

    plotTitle <- main

    ## for MIRROR functionality
    arg.list <- formals(xpose.plot.histogram)
    arg.names <- names(arg.list)
    new.arg.list <- vector("list",length(arg.names))
    names(new.arg.list) <- arg.names
    for (argnam in arg.names){
      if (argnam=="..."){
        next
      }
      tmp <- get(argnam)
      if (is.null(tmp)){
      } else {
        new.arg.list[[argnam]]=tmp
      }
    }
    if (mirror){
      create.mirror(xpose.plot.histogram,
                         new.arg.list,mirror,plotTitle,...)
    } else { # end if mirror

      
      ## x-label
      if(!is.null(x)) {
        if(length(x)> 1) {
          xlb <- NULL
        } else {
          if(is.null(xlb)) {
            xlb <- xlabel(x,object)
          }
        }
      }
      
      ##Get data
      if(!is.null(samp)) {
        data <- SData(object,inclZeroWRES,onlyfirst=onlyfirst,
                      subset=subset,samp=samp)
      } else {
        data <- Data(object,inclZeroWRES,onlyfirst=onlyfirst,subset=subset)
      }
      
      ## Strip "missing" data
      ##data <- subset(data, get(x) != object@Prefs@Miss)
      
      ##if(any(is.null(data))) return("The subset expression is invalid!\n")

      ## Make sure by is a factor if requested
      if(!is.null(by) && force.by.factor) {
        for(b in by) {
          data[,b] <- as.factor(data[,b])
        }
      }

      ## Check to see if more that one x-variable
      if(length(x) > 1) {
        reps <-c(xvardef("id",object),xvardef("idlab",object),
                 xvardef("wres",object))
        
        if(!is.null(by)) reps <- c(reps,by)
        data <- xpose.stack(data,object,x,reps)
        
        object <- new("xpose.data",
                      Runno=object@Runno,
                      Data = NULL)
        Data(object) <- data

        onlyfirst = FALSE
        if(is.null(by)) {
          by <- "ind"
        } else {
          by <- c("ind",by)
        }

        x <- "values"
        scales=list(relation="free")
      }

      ## Strip "missing" data
      data <- subset(data, get(x) != object@Prefs@Miss)
      
      if(any(is.null(data))) return("The subset expression is invalid!\n")

      ## Collect the basic plot formula
      bb <- NULL
      if(any(is.null(by))) { ## No conditioning
        formel <- paste("~",x,sep="")
      } else {
        for(b in by) {
          bb <- c(bb,xlabel(b,object))
          
          if(!is.factor(data[,b])) {
            data[,b] <- equal.count(data[,b],number=shingnum,overl=shingol)
          } else {
            
            if(any(!is.null(ordby))) {
              data[,b] <- reorder(data[,b],data[,ordby],byordfun)
            }
            
            if(names(data[,b,drop=F])!="ind") {
              levels(data[,b]) <-
                paste(xlabel(names(data[,b,drop=F]),object),":",   ## Needs to be fixed
                      levels(data[,b]),sep="")
            }
            
          }
        }
        bys    <- paste(by,collapse="*")
        formel <-  paste("~",x,"|",bys,sep="")
      }

      if(missing(strip)) {
        strip <- function(var.name,...)
          strip.default(var.name=bb,strip.names=c(F,T),...)
      }

      
      
      xvarnam <- x

      if(!is.null(x.cex)) {
        if (is.list(xlb)){
          xlb$cex=x.cex
        } else {
          xlb <- list(xlb,cex=x.cex)
        }
      }
      if(!is.null(y.cex)) {
        if (is.list(ylb)){
          ylb$cex=y.cex
        } else {
          ylb <- list(ylb,cex=y.cex)
        }
      }
      
      if(is.null(main)) {
      } else {
        if(!is.null(main.cex)) {
          if (is.list(main)){
            main$cex=main.cex
          } else {
            main <- list(main,cex=main.cex)
          }
        }
      }

      if(missing("type")) {
        if (length(levs <- unique(data[,x])) <= object@Prefs@Cat.levels) {
          type <- "count"
          ylb <- "Count"
        }
      }


      xplot <- histogram(formula(formel),data,
                         ...,
                         obj=object,
                         prepanel = function(x,bins.per.panel.equal = TRUE,...) {
                           if(length(levs <- unique(x)) <= object@Prefs@Cat.levels) {
                             xlim <- as.character(sort(levs))
                             return(list(xlim=xlim))
                           } else {
                             xlim <- range(x)
                             
                             if(!bins.per.panel.equal){
                                 nint      <- round(log2(length(x))+1)
                                 endpoints <- range(x[!is.na(x)])
                                 breaks <- do.breaks(endpoints, nint)
                                 hdat <- hist(x, breaks=breaks, plot=F)
                                 ddat <- density(x,na.rm=T)
                                 ylim <- c(0,  max(hdat$density,ddat$y))
                                 return(list(xlim=xlim,ylim=ylim))
                             }
                                 return(list(xlim=xlim))
                           }
                         },
                         panel=xpose.panel.histogram,
                         aspect=aspect,
                         ylab=ylb,
                         xlab=xlb,
                         type=type,
                         scales=scales,
                         main=main,
                         xvarnam=xvarnam,
                         hidlty = hidlty,
                         hidcol = hidcol,
                         hidlwd = hidlwd,
                         hiborder = hiborder,
                         hilty = hilty,
                         hilwd = hilwd,
                         strip=strip,
                         hicol = hicol)
      
      
      return(xplot)
    }
  }
