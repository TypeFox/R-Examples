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

"xpose.plot.default" <-
  function(x,y,object,
           inclZeroWRES = FALSE,
           onlyfirst    = FALSE,
           samp         = NULL,
           panel        = xpose.panel.default,
           groups       = object@Prefs@Xvardef$id,
           ids          = object@Prefs@Graph.prefs$ids,
           logy         = FALSE,
           logx         = FALSE,
           yscale.components= "default",#function(...) yscale.components.default(...),
           xscale.components= "default",#function(...) xscale.components.default(...),
           
           aspect       = object@Prefs@Graph.prefs$aspect,
           funx         = NULL,
           funy         = NULL,
           iplot        = NULL,
           
           ## Prediction interval settings
           PI           = NULL,
           
           ## Conditioning settings
           by           = object@Prefs@Graph.prefs$condvar,
           force.by.factor = FALSE,
           ordby        = object@Prefs@Graph.prefs$ordby,
           byordfun     = object@Prefs@Graph.prefs$byordfun,
           shingnum     = object@Prefs@Graph.prefs$shingnum,
           shingol      = object@Prefs@Graph.prefs$shingol,
           by.interval  = NULL,
           ##par.strip.text=trellis.par.get("add.text"),
           ##mirror.par.strip.text=trellis.par.get("add.text"),
           strip = function(...){
             strip.default(...,strip.names=c(TRUE,TRUE))
           },
           use.xpose.factor.strip.names=TRUE,
           ##strip.nams=T,
           ##strip=strip.custom(strip.names=c(T,T)),
           ##par.strip.text = mirror.par.strip.text),
           ##par.strip.text = trellis.par.get("add.text"),
           ##par.strip.text=NULL,
           
           ## Subset stuff
           subset       = xsubset(object),
           
           autocorr=FALSE,
           
           ## Axes and titles
           main         = xpose.create.title(x,y,object,subset,funx,funy,...),
           #main         = NULL,
           xlb          = xpose.create.label(x,object,funx,logx,autocorr.x=autocorr,...),
           ylb          = xpose.create.label(y,object,funy,logy,autocorr.y=autocorr,...),
           ##xlb          = ifelse((length(x)>1),"Value",xlabel(x,object)),
           ##ylb          = ifelse((length(y)>1),"Value",xlabel(y,object)),
           scales       = list(),           
           
           ## Superpose smooth
           suline       = object@Prefs@Graph.prefs$suline,
           
           ## Categorical stuff
           bwhoriz      = object@Prefs@Graph.prefs$bwhoriz,
           
           ## Dilution stuff
           dilution     = FALSE,
           dilfrac      = object@Prefs@Graph.prefs$dilfrac,
           diltype      = object@Prefs@Graph.prefs$diltype,
           dilci        = object@Prefs@Graph.prefs$dilci,
           seed         = NULL,
           
           
           
           mirror       = FALSE,
           max.plots.per.page=4,
           mirror.aspect="fill",
           pass.plot.list=FALSE,
           x.cex        = NULL,
           y.cex        = NULL,
           main.cex     = NULL,
           mirror.internal=list(strip.missing=missing(strip)),
           ...
  ) {
    
    ## CHecks if use.xpose.factor.strip.names is a logical a length 1
    if (!(class(use.xpose.factor.strip.names)=="logical" & 
            length(use.xpose.factor.strip.names)==1)){
      stop("The provided use.xpose.factor.strip.names argument is not a logical of length 1")
    }
    
    plotTitle <- main
    
    ## for MIRROR functionality
    arg.list <- formals(xpose.plot.default)
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
      if(is.null(object@Nsim)) {
        cat(paste("The current Xpose database does not have any simulation data.\n"))
        cat(paste("The mirror option cannot be used.\n"))
        return(NULL)
      }
      create.mirror(xpose.plot.default,
                    new.arg.list,mirror,plotTitle,...)
    } else { # end if mirror
      
      
      ##Get data
      if(any(is.null(iplot))) {
        if(!is.null(samp)) {
          #cat(samp)
          data <- SData(object,inclZeroWRES,onlyfirst=onlyfirst,
                        subset=subset,samp=samp)
        } else {
          data <- Data(object,inclZeroWRES,onlyfirst=onlyfirst,subset=subset)
        }
      } else {
        data <- Data(object,inclZeroWRES,onlyfirst=onlyfirst,subset=NULL)
      }
      
      ## Strip "missing" data
      data <- subset(data, get(x) != object@Prefs@Miss)
      data <- subset(data, get(y) != object@Prefs@Miss)
      
      if(any(is.null(data))) return("The subset expression is invalid.")
      
      ## Make sure by is a factor if requested
      if(!is.null(by) && force.by.factor) {
        for(b in by) {
          data[,b] <- as.factor(data[,b])
        }
      }
      
      ## Sort out dilution
      dilsubset <- TRUE
      dilname   <- NULL
      if(dilution) {
        if(is.null(diltype)) { # Standard random dilution
          data <- create.rand(data,object,dilfrac,seed=seed)
          if(is.null(seed)) {
            dilsubset <- parse(text="Rnoseed==0")
            dilname   <- "Rnoseed"
          } else {
            dilsubset <- parse(text=paste("R",seed,"==0",sep=""))
            dilname   <- paste("R",seed,"==0",sep="")
          }
        } else {               # Stratified random dilution
          data <-create.strat.rand(data,object,x,y,dilfrac,dilci,seed=seed)
          if(is.null(seed)) {
            dilsubset <- parse(text="RSnoseed==0")
            dilname   <- "RSnoseed"
          } else {
            dilsubset <- parse(text=paste("RS",seed,,"==0",sep=""))
            dilname   <- paste("RS",seed,,"==0",sep="")
          }
        }
      }
      
      ## Check to see if x and y are both longer than 1
      if(length(x)>1 && length(y)>1) {
        cat("x and y can not both be longer than 1\n")
        return()
      }
      
      
      ## Check to see if more than one x-variable
      if(length(x) > 1) {
        reps <-c(xvardef("id",object),xvardef("idlab",object),
                 xvardef("wres",object),y,groups)
        if(!is.null(dilname)) reps <- c(reps,dilname)
        
        if(!is.null(by)) reps <- c(reps,by)
        #data <- stack.xpose(data,object,x,reps)
        data <- xpose.stack(data,object,x,reps)
        object <- new("xpose.data",
                      Runno=object@Runno,
                      Data = NULL,
                      Prefs = object@Prefs)
        
        Data(object) <- data
        #cat(object@Prefs@Graph.prefs$type)
        
        if(is.null(main.cex)) main.cex <- 0.9
        onlyfirst = FALSE
        if(is.null(by)) {
          by <- "ind"
        } else {
          by <- c("ind",by)
        }
        
        x <- "values"
        
        ## If scales is longer than one then the users has supplied it
        ##as an argument.
        if(length(scales)==0) {
          scales=list(x=list(relation="free"))
        }
      }
      
      ## Check to see if more than one y-variable
      if(length(y) > 1) {
        reps <- c(object@Prefs@Xvardef["id"],
                  object@Prefs@Xvardef["idlab"],
                  xvardef("wres",object),x,groups)
        if(!is.null(dilname)) reps <- c(reps,dilname)
        
        if(!is.null(by)) reps <- c(reps,by)
        #data <- stack.xpose(data,object,y,reps)
        data <- xpose.stack(data,object,y,reps)
        object <- new("xpose.data",
                      Runno=object@Runno,
                      Data = NULL,
                      Prefs = object@Prefs)
        
        Data(object) <- data
        
        if(is.null(main.cex)) main.cex <- 0.9
        onlyfirst = FALSE
        
        if(is.null(by)) {
          by <- "ind"
        } else {
          by <- c("ind",by)
        }
        
        y <- "values"
        
        ## If scales is longer than one then the users has supplied it
        ##as an argument.
        if(length(scales)==0) {
          scales=list(y=list(relation="free"))
        }
      }
      
      
      ## Collect the basic plot formula
      bb <- NULL
      groups <- groups
      if(any(is.null(by))) {
        if(bwhoriz) {
          formel <- paste(x,"~",y,sep="")
        } else {
          formel <- paste(y,"~",x,sep="")
        }
      } else {
        for(b in by) {
          
          ##b <- by[bs]
          
          bb <- c(bb,xlabel(b,object))
          
          if(!is.factor(data[,b])) {
            if(is.null(by.interval)){
              data[,b] <- equal.count(data[,b],number=shingnum,overl=shingol)
            } else {
              data[,b] <- shingle(data[,b],intervals=by.interval)
            }
          } else {
            
            if(any(!is.null(ordby))) {
              data[,b] <- reorder(data[,b],data[,ordby],byordfun)
            }
            
            if(names(data[,b,drop=F])!="ind") {
              if(use.xpose.factor.strip.names){
                levels(data[,b]) <-
                  paste(xlabel(names(data[,b,drop=F]),object),":",   ## Needs to be fixed
                        levels(data[,b]),sep="")
              }
            }
          }
        }
        bys    <- paste(by,collapse="*")
        if(bwhoriz) {
          formel <-  paste(x,"~",y,"|",bys,sep="")
        } else {
          formel <-  paste(y,"~",x,"|",bys,sep="")
        }
      }
      
      if(missing(strip)) {
        strip <- function(var.name,...)
          strip.default(var.name=bb,strip.names=c(F,T),...)
      }
      
      ## Check to see if panel.superpose should be used
      if(any(!is.null(groups))) groups <- data[,groups]
      
      ## CHeck to see if a superpose smooth is to be used.
      if(!is.null(suline)) {
        suline <- data[,suline]
      }
      
      ## Check for id-numbers as plotting symbols
      ##if(!is.null(ids)) ids <- data[,xvardef("idlab",object)]
      if(ids){
        ids <- data[,xvardef("idlab",object)]
      } else {
        ids <- NULL
      }
      
      ## Apply function to x-variable
      if(!is.null(funx)) {
        data[,x] <- do.call(funx,list(data[,x]))
      }
      
      ## Apply function to y-variable
      if(!is.null(funy)) {
        data[,y] <- do.call(funy,list(data[,y]))
        ##         if(!is.null(ylb[1])){
        ##           ##if(ylb[1]==xlabel(y,object)) {
        ##           if(missing(ylb)) {
        ##             if (fun=="abs"){
        ##               ylb <- paste("|",ylb,"|",sep="")
        ##             } else {
        ##               ylb <- paste(fun,"(",ylb,")",sep="")
        ##             }
        ##           }
        ##         }
      }
      
      ## Sort out the scales
      yscale.components.defined <- T
      xscale.components.defined <- T
      
      if(!is.function(yscale.components)){
        if(!is.na(match(yscale.components,"default"))) {
          yscale.components= function(...) yscale.components.default(...)
          yscale.components.defined <- F
        }
      }
      
      if(!is.function(xscale.components)){
        if(!is.na(match(xscale.components,"default"))) {
          xscale.components= function(...) xscale.components.default(...)
          xscale.components.defined <- F
        }
      }
      
      if(logy) {
        scales$y$log <- TRUE
        if(!yscale.components.defined){
          yscale.components=xpose.yscale.components.log10
        }
      }
      if(logx) {
        scales$x$log <- TRUE
        if(!xscale.components.defined){
          xscale.components=xpose.xscale.components.log10
        }
      }
      
      
      xvarnam <- x
      yvarnam <- y
      
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
      
      
      ## for autocorrelation (not working completely yet)
      if(autocorr){
        auto.ids <- unique(data[[xvardef("id",object)]])
        auto.n <- 0
        xplt1 <- 0
        xplt2 <- 0
        xgrps <- 0
        for(i in 1:length(auto.ids)) {
          i <- 1
          seli <- data[[xvardef("id",object)]] == ids[i]
          nobs <- length(data[[x]][seli])
          xplt <- matrix(data[[x]][seli], 1, nobs)
          if(nobs > 1) {
            for(j in 1:(nobs - 1)) {
              auto.n <- auto.n + 1
              xplt1[auto.n] <- xplt[1, j]
              xplt2[auto.n] <- xplt[1, j + 1]
              xgrps[auto.n] <- auto.ids[i]
            }
          }
        }
        
        #xlb <- paste(xlb,"(i)",sep="")
        #ylb <- paste(ylb,"(i+1)",sep="")
        
        #x <- xplt1
        #y <- xplt2
        #groups <- xgrps
      }
      
      xplot <- xyplot(formula(formel),data,obj=object,
                      prepanel = function(x,y) {
                        xlim <- NULL
                        ylim <- NULL
                        ret <- list()
                        if(is.factor(x)){#length(levs <- unique(x)) < object@Prefs@Cat.levels) {
                          if(length(grep("[A-Z,a-z]",levels(x)))==0) {
                            xlim <- as.character(sort(as.numeric(levels(x))))
                          } else {
                            xlim <- sort(levels(x))
                          }
                        } else {
                          #xlim <- range(x)
                        }
                        ret[["xlim"]] <- xlim
                        if(is.factor(y)){#length(levs <- unique(x)) < object@Prefs@Cat.levels) {
                          if(length(grep("[A-Z,a-z]",levels(y)))==0) {
                            ylim <- as.character(sort(as.numeric(levels(y))))
                          } else {
                            ylim <- sort(levels(y))
                          }
                        } else {
                          #ylim <- range(y)
                        }
                        ret[["ylim"]] <- ylim
                        #list(xlim=xlim,ylim=ylim)
                        return(ret)
                      },
                      onlyfirst = onlyfirst,
                      samp   = samp,
                      panel = panel,
                      strip = strip,
                      ##par.strip.text = par.strip.text,
                      groups=groups,
                      inclZeroWRES=inclZeroWRES,
                      PI    = PI,
                      logy=logy,
                      logx=logx,
                      xscale.components=xscale.components,
                      yscale.components=yscale.components,
                      xvarnam = xvarnam,
                      yvarnam = yvarnam,
                      ids   = ids,
                      main=main,
                      xlab=xlb,
                      ylab=ylb,
                      aspect=aspect,
                      suline=suline,
                      bwhoriz = bwhoriz,
                      subset=eval(dilsubset),
                      scales=scales,
                      iplot=iplot,
                      autocorr=autocorr,
                      #autocorr=FALSE,
                      PI.subset=subset,
                      #drop.unused.levels=FALSE,
                      ...)
      return(xplot)
    }
  }

