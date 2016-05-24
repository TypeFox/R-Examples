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

## Added by Andrew Hooker
## 21/12/2005

create.mirror <- function(fun,arg.list,mirror,plotTitle,
                               fix.y.limits=TRUE,
                               fix.x.limits=TRUE,
                               ...){
  
  ## what sort of mirror do we have?
  if(mirror) {      
    if(!is.logical(mirror)) {
      if(mirror != 1 && mirror !=3) {
        cat("The mirror should either be logical, 1 or 3!\n")
        invisible()
        return(NULL)
      }
    } else {
      mirror <- 1
    }
  }

  arg.list$mirror=FALSE
  arg.list$...=NULL
  
  if (length(grep("par.strip.text",deparse(match.call())))>0){
    par.strip.text.exists <- TRUE
  } else {
    par.strip.text.exists <- FALSE
  }

  if (length(grep("auto.key",deparse(match.call())))>0){
    auto.key.exists <- TRUE
  } else {
    auto.key.exists <- FALSE
  }


  ## if strip was not supplied by user and exists then get rid of argument
  if(!is.null(arg.list$strip)){
    if(!is.null(arg.list$mirror.internal$strip.missing)){
      if(arg.list$mirror.internal$strip.missing){
        arg.list$strip=NULL
      }
    }
  }
  
  ## Set the seed number and decide what simulated data sets to plot
  if(!is.null(arg.list$seed)) set.seed(arg.list$seed)
  rand.samp <- sample(1:arg.list$object@Nsim,mirror)

  ## size of labels
  if(mirror==3) {
    if(is.null(arg.list$y.cex)) arg.list$y.cex=0.6
    if(is.null(arg.list$x.cex)) arg.list$x.cex=0.6
  }          
  if(mirror==1) {
    arg.list$y.cex=0.8
    arg.list$x.cex=0.8
  }

  if(mirror==1){
    if(is.null(arg.list$scales$cex)) arg.list$scales$cex=0.7
    if(is.null(arg.list$scales$tck)) arg.list$scales$tck=0.7
    #if(is.null(arg.list$scales$cex)) arg.list$scales=list(cex=0.7,tck=0.7)
    if(is.null(arg.list$main.cex)) arg.list$main.cex =0.9
    if(!par.strip.text.exists) arg.list$par.strip.text$cex =0.9
    if(!is.null(arg.list$auto.key)){
      if(arg.list$auto.key=="Default") arg.list$auto.key <- list(cex=0.8)
    }
  }
  
  if(mirror==3){
    if(is.null(arg.list$scales$cex)) arg.list$scales$cex=0.5
    if(is.null(arg.list$scales$tck)) arg.list$scales$tck=0.5
    #if(is.null(arg.list$scales$cex)) arg.list$scales=list(cex=0.5,tck=0.5)
    if(is.null(arg.list$main.cex)) arg.list$main.cex =0.7
    if(!par.strip.text.exists) arg.list$par.strip.text$cex =0.6
    if(!is.null(arg.list$auto.key)){
      if(arg.list$auto.key=="Default") arg.list$auto.key <- list(cex=0.6)
    }
  }

  ## aspect of graphs
  arg.list$aspect=arg.list$mirror.aspect

  
  total.plot.number <- mirror+1
  full.plot.list <- vector("list",total.plot.number)

  xlb <- arg.list$xlb
  for (j in 0:mirror){
    
    ## The plot titles
    if (j==0) arg.list$main  <- "Original data"
    if (j!=0) arg.list$main <- paste("Simulated data (#",rand.samp[j],")",sep="")

    ## The axis labels
    arg.list$xlb  <- xlb
    if (mirror==3 & (j==1 | j==2) ) arg.list$xlb  <- " "
    if (mirror==1 & j==0) arg.list$xlb  <- " "

    
    ## the simulation samples
    if (j==0) arg.list$samp=NULL
    if (j!=0) arg.list$samp=rand.samp[j]
    

    arg.string=NULL
    for(argnam in names(arg.list)){
      if (is.null(arg.string)){
        arg.string="list("
      }
      arg.string=paste(arg.string,argnam,"=arg.list$",argnam,",",sep="")
    }
    arg.string=paste(arg.string,"...)",sep="")    
    evaluated.arg.string <- eval(parse(text=arg.string))
    plot.val <- do.call(fun,evaluated.arg.string)

    full.plot.list[[j+1]] <- plot.val
  }

  
  if(fix.x.limits) xlimits=c()
  if(fix.y.limits) ylimits=c()
  
 
  
  
  for (j in 0:mirror){
    xplot <- full.plot.list[[j+1]]
    if(fix.x.limits) xlimits <- c(xlimits,xplot$x.limits)
    if(fix.y.limits) ylimits <- c(ylimits,xplot$y.limits)
  }
  if(!is.numeric(xlimits)) fix.x.limits <- FALSE
  if(!is.numeric(ylimits)) fix.y.limits <- FALSE
  if(fix.x.limits) new.xlimits <- c(min(xlimits),max(xlimits))
  if(fix.y.limits) new.ylimits <- c(min(ylimits),max(ylimits))
  
  for (j in 0:mirror){
    if(fix.x.limits) full.plot.list[[j+1]]$x.limits <- new.xlimits
    if(fix.y.limits) full.plot.list[[j+1]]$y.limits <- new.ylimits
  }
  
  ##if(is.null(arg.list$max.plots.per.page)) arg.list$max.plots.per.page=4
  if(arg.list$pass.plot.list){
    return(full.plot.list)
  } else {

    obj <- xpose.multiple.plot(full.plot.list,plotTitle=plotTitle,
                               #prompt=FALSE,#arg.list$prompt,
                               #max.plots.per.page=arg.list$max.plots.per.page,
                               mirror=mirror,
                               ...)
    return(obj)
  }
}
