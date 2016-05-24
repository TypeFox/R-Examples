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

"xpose.plot.splom" <-
  function(plist, object,
           varnames=NULL,
           main = "Scatterplot Matrix",
           xlb = NULL,
           ylb = NULL,
           scales = list(),
           onlyfirst=TRUE,
           inclZeroWRES=FALSE,
           subset = xsubset(object),
           by           = object@Prefs@Graph.prefs$condvar,
           force.by.factor=FALSE,
           include.cat.vars = FALSE,
           ordby     = NULL,
           byordfun  = object@Prefs@Graph.prefs$byordfun,
           shingnum  = object@Prefs@Graph.prefs$shingnum,
           shingol   = object@Prefs@Graph.prefs$shingol,
           strip = function(...)
           strip.default(...,strip.names=c(TRUE,TRUE)),
           #par.strip.text=trellis.par.get("add.text"),
           groups = NULL,
           ids = object@Prefs@Graph.prefs$ids,
           smooth       = TRUE,
           lmline = NULL,
           panel        = xpose.panel.splom,
           aspect = object@Prefs@Graph.prefs$aspect,
                                        #varname.cex=NULL,
                                        #axis.text.cex=NULL,
           ## mirror stuff
           samp=NULL,
           max.plots.per.page=4,
           mirror       = FALSE,
           mirror.aspect="fill",
           pass.plot.list=FALSE,
           x.cex=NULL,
           y.cex=NULL,
           main.cex=NULL,
           mirror.internal=list(strip.missing=missing(strip)),
           ...) {

    
    plotTitle <- main

    ## for MIRROR functionality
    arg.list <- formals(xpose.plot.splom)
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
      create.mirror(xpose.plot.splom,
                         new.arg.list,mirror,plotTitle,...)
    } else { # end if mirror
      
      ##Get data
                                        #data <- object@Data[, xvardef("parms", object), drop = F]
                                        #mlist <- c(plist, by, groups)
                                        #data <- Data(object,inclZeroWRES,onlyfirst=onlyfirst,subset=subset)[, mlist, drop = F]
      if(!is.null(samp)) {
        data <- SData(object,inclZeroWRES=inclZeroWRES,onlyfirst=onlyfirst,
                      subset=subset,samp=samp)
      } else {
        data <- Data(object,inclZeroWRES=inclZeroWRES,onlyfirst=onlyfirst,subset=subset)
      }

      ## Strip "missing" data
      for (i in plist) {
        data <- subset(data, get(i) != object@Prefs@Miss)
      }

      if(any(is.null(data))) return("The data or subset expression is invalid.")
      
      ## if the parameter or variable in the list has only one value don't plot it
      remove.from.plist=c()

      for (i in 1:length(plist)) {
        if(!is.factor(data[,plist[i]])){
          if(length(unique(data[,plist[i]])) < 2){
            remove.from.plist=c(remove.from.plist,i)
            cat(paste(plist[i],
                      "has only one value and will not be\n",
                      "shown in the scatterplot\n"))
          } 
        } else {
          if(!include.cat.vars){
            remove.from.plist=c(remove.from.plist,i)
            cat(paste(plist[i],
                      "is categorical and will not be\n",
                      "shown in the scatterplot\n"))
          } else {
            if(length(levels(data[,plist[i]])) < 2){
              remove.from.plist=c(remove.from.plist,i)
              cat(paste(plist[i],
                        " has only one value and will not be\n",
                        "shown in the scatterplot\n",sep=""))
            }
          }
        }
      }
      if(length(remove.from.plist)>0){
        plist <- plist[-remove.from.plist]
        if(!is.null(varnames))  varnames <- varnames[-remove.from.plist]
      }
      
      if(is.null(varnames)) {
        varnames <- c()
        for (i in plist) {
          varnames <- c(varnames, xlabel(i, object))
        }
      }
      
      ## Make sure by is a factor if requested
      if(!is.null(by) && force.by.factor) {
        for(b in by) {
          data[,b] <- as.factor(data[,b])
        }
      }
      
      ## Collect the basic plot formula
      bb <- NULL

      if(any(is.null(by))) {
        formel <- paste("~data[, plist]", sep="") 
      } else {
        for(b in by) {
                                        #b <- by[bs]

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
        formel <-  paste("~data[, plist] | ", bys, sep="")
      }

      if(missing(strip)) {
        strip <- function(var.name,...)
          strip.default(var.name=bb,strip.names=c(F,T),...)
      }

      ## Check to see if panel.superpose should be used
      if(any(!is.null(groups))) groups <- data[,groups]

      ## CHeck to see if a superpose smooth is to be used.
      suline <- NULL
      if(!is.null(suline)) {
        suline <- data[,suline]
      }
      
      ## Check for id-numbers as plotting symbols
      if(!is.null(ids)) ids <- data[,xvardef("idlab",object)]

                                        #cat(formel)
                                        #readline()
                                        #browser()

      ##     if(length(plist)>7) {
      ##       varname.cex=0.6
      ##       axis.text.cex=0.6
      ##     }

      
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

      xplot <- splom(formula(formel), data, obj=object,
                                        #prepanel.limits = function(x) 
                                        #  if (is.factor(x)) levels(x) else
                                        #    extend.limits(range(as.numeric(x), na.rm = TRUE)),
                     varnames=varnames,
                     onlyfirst = onlyfirst,
                     panel=panel,
                     strip = strip,
                                        #par.strip.text = par.strip.text,
                     groups=groups,
                     inclZeroWRES=inclZeroWRES,
                     ids   = ids,
                     main=main,
                                        #aspect=aspect,
                     smooth=smooth,
                     lmline = lmline,
                     ylab = ylb,
                     xlab = xlb,
                     scales = scales,
                                        #varname.cex=varname.cex,
                                        #axis.text.cex=axis.text.cex,
                     ...)
      return(xplot)
      
    }

  }
