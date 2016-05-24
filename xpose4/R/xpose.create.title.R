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

xpose.create.title <- function(x,y,object,subset=NULL,funx=NULL,funy=NULL,
                               no.runno=FALSE,...){

  vs.label <- " vs. "
  
  x.name <- xlabel(x,object)
  if(length(x) > 1) {
    x.name <- NULL
    for(xx in x){
      if (is.null(x.name)){
        x.name <- paste(xlabel(xx,object))
      } else {
        x.name <- paste(x.name,"/",xlabel(xx,object))
      }
    }
    vs.label <- paste(vs.label,"\n")
    ##x.name <- paste(x.name,"\n")
  } # end length>1
  
  
  y.name <- xlabel(y,object)
  if(length(y) > 1) {
    y.name <- NULL
    for(yy in y){
      if (is.null(y.name)){
        y.name <- paste(xlabel(yy,object))
      } else {
        y.name <- paste(y.name,"/",xlabel(yy,object))
      }
    }
    vs.label <- paste("\n",vs.label)
  }# end length>1
  
  
  main.runno <- ifelse(!no.runno,
                       paste(" (Run ",object@Runno,")",sep=""),
                       "")


  if(is.null(x.name)) vs.label <- NULL
  if(!is.null(funx)) {
    if(is.null(x.name)) {
      main.xname <- x.name
    } else {
      main.xname <- paste(funx,"(",x.name,")",sep="")
      if (funx=="abs"){main.xname <- paste("|",x.name,"|",sep="")}
    }
  } else {
    main.xname <- x.name
  }
  
  if(!is.null(funy)) {
    main.yname <- paste(funy,"(",y.name,")",sep="")
    if (funy=="abs"){main.yname <- paste("|",y.name,"|",sep="")}
  } else {
    main.yname <- y.name
  }

  main <- paste(main.yname, vs.label, main.xname, main.runno, sep="")
  
  if (!is.null(subset)){
    main <- paste(main,"\n[",subset,"]",sep="")
  }
  return(main)
}
 
xpose.create.title.text <- function(x,y,text,object,subset,text2=NULL,...){
  main <- xpose.create.title(x,y,object,subset=subset,...)
  main <- paste(text,main,text2)
  return(main)
}

xpose.create.title.hist <- function(x,object,subset,...){
  main <- paste("Distribution of ",xlabel(x,object),
                " (Run ",object@Runno,")",sep="")
  if (!is.null(subset)){
    main <- paste(main,"\n[",subset,"]",sep="")
  }
  return(main)
}

xpose.create.label <- function(x,object,fun,logx,
                               autocorr.x=FALSE,
                               autocorr.y=FALSE,...){

  x.label <- ifelse((length(x)>1),"Value",xlabel(x,object))

  tot.x.label <- x.label
  #if(logx) tot.x.label <- paste("log(",tot.x.label,")",sep="")
  if(!is.null(fun)) {
    tot.x.label <- paste(fun,"(",x.label,")",sep="")
    if (fun=="abs"){
      tot.x.label <- paste("|",x.label,"|",sep="")
    }
  }
  if(autocorr.x) tot.x.label <- paste(tot.x.label,"(i)",sep="")
  if(autocorr.y) tot.x.label <- paste(tot.x.label,"(i+1)",sep="")
  
  return(tot.x.label)
}

xpose.multiple.plot.title <-
  function(object,
           plot.text,
           subset=xsubset(object),
           main="Default",
           no.runno=FALSE,
           ...){
    if (is.null(main)){
      plotTitle <- NULL
    } else {
      if(!is.na(match(main,"Default"))) {
        plotTitle <- paste(plot.text," (Run ",object@Runno, ")", sep="")
        if (no.runno) plotTitle <- paste(plot.text, sep="")
        if (!is.null(subset)){
          plotTitle <- paste(plotTitle,"\n[",subset,"]",sep="")
        }
      } else {
        plotTitle <- main
      }
    }
    return(plotTitle)
  }
