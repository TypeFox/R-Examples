# Jasmin Straube, Queensland Facility of Advanced Bioinformatics
# Part of this script was borrowed the graphics package.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU Moleculesral Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Moleculesral Public License for more details.
#
# You should have received a copy of the GNU Moleculesral Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#' Plot of \code{lmmspline} object
#' 
#' Plots the raw data, the mean and the fitted or derivative information of the \code{lmmspline} object.
#' 
#' @import ggplot2
#' @importFrom stats spline na.omit
#' @param x An object of class \code{lmmspline}.
#' @param y \code{character} or \code{numeric} value. Determining which feature should be plotted can be either the index or the name of the feature. 
#' @param data alternative \code{matrix} or \code{data.frame} containing the original data for visualisation purposes.
#' @param time alternative \code{numeric} indicating the sample time point. Vector of same length as row length of data for visualisation purposes.
#' @param smooth an optional \code{logical} value. Default \code{FALSE}, if \code{TRUE} smooth representation of the fitted values. 
#' @param mean alternative \code{logical} if the mean should be displayed. By default set to \code{TRUE}.
#' @param \ldots Additional arguments which are passed to \code{plot}.
#' @return plot showing raw data, mean profile and fitted profile. 
#' @examples 
#' \dontrun{
#' data(kidneySimTimeGroup)
#' # running for samples in group 1
#' G1 <- which(kidneySimTimeGroup$group=="G1")
#' testLmmspline <- lmmSpline(data=kidneySimTimeGroup$data[G1,],
#'                  time=kidneySimTimeGroup$time[G1],
#'                  sampleID=kidneySimTimeGroup$sampleID[G1],keepModels=T)
#'                  
#'          
#' plot(testLmmspline, y=2)
#' plot(testLmmspline, y=2, smooth=TRUE)
#' # Don't keep the models to improve memory usage 
#' testLmmspline <- lmmSpline(data=kidneySimTimeGroup$data[G1,],
#'                  time=kidneySimTimeGroup$time[G1],
#'                  sampleID=kidneySimTimeGroup$sampleID[G1],keepModels=F)
#'                  
#' #plot only the fitted values
#' plot(testLmmspline, y=2)
#' #plot fitted values with original data
#' plot(testLmmspline, y=2, data=kidneySimTimeGroup$data[G1,], time=kidneySimTimeGroup$time[G1])
#' 
#' }

#' @method plot lmmspline
#' @export
plot.lmmspline <- function(x, y, data,time, smooth,mean,...){
  
  if(length(y)>1)
    stop('Can just plot a single feature.')
  name <- ""
  if(sum(is.na(x@predSpline[y,]))==length(x@predSpline[y,]))
    stop('Error plotting molecule')
  if(missing(smooth))
    smooth <- F
  if(class(y)=='numeric'){
    name <- rownames(x@predSpline)[y]
  }
  if(class(y)=='character'){
    nam <- rownames(x@predSpline)
    if(sum(nam%in%y)>0){
      name <- y
      y <- which(nam%in%y)
     
   
    }else{
      stop(paste('Could not find feature',y,'in rownames(x@pred.spline).'))
    }
  }
  if(length(x@models)>0){
    model <- x@models[[y]]
  if(x@derivative){
    p <- plotLmmsDeriv(model, smooth=smooth,data=x@predSpline[y,],name,...) 
  }else{
     p <- plotLmms(model,smooth=smooth,name,mean,...)
  }
  }else{
    t <- as.numeric(colnames(x@predSpline))
    
    if(x@derivative){
      p <- plotdataLMMS(t,x@predSpline[y,],smooth,name,data,time,y,mean,...)
    }else{

      p <- plotdataLMMS(t,x@predSpline[y,],smooth,name,data,time,y,mean,...)
  }
}
return(p)
}

plotdataLMMS <- function(x,y,smooth,name,data,time,mol,mean,...){
  Time <- Intensity <- Model <- NULL
  if(sum(is.na(x))==length(x))
    return('Error plotting molecule')
  g <- ggplot()+ggtitle(name)
  if(missing(mean))
    mean <- FALSE
  if(missing(data)|missing(time)){  
    if(smooth){
      spl <- spline(x = x, y = unlist(y), n = 500, method = "natural")
      s <- data.frame(Time=spl$x,Intensity=spl$y,Model="Smooth")
    }else{    
      s <- data.frame(Time = x,Intensity = unlist(y),Model="Fitted")
    }
    g <- g+ geom_line(aes(x = Time,y=Intensity,linetype=Model),data = s,size=1)
    if(mean)
      g <- g + stat_summary(aes(x=Time,y=Intensity,linetype='Mean'),size=1,data=s,fun.y=function(x)mean(x), geom="line",na.rm = T)
    
  }else{
    s <- data.frame(Time = time,Intensity = data[,mol],Model="Mean")
    g <- g+ geom_point(aes(x = Time,y=Intensity),alpha=0.5,data = s,size=3,na.rm = T)
    if(mean)
      g <- g + stat_summary(aes(x=Time,y=Intensity,linetype=Model),size=1,data=s,fun.y=function(x)mean(x), geom="line",na.rm = T)
    if(smooth){      
      spl <- spline(x = x, y = unlist(y), n = 500, method = "natural")
      
      s <- data.frame(Time = spl$x,Intensity = spl$y,Model="Smooth")
    }else{
      s <- data.frame(Time = x,Intensity=unlist(y),Model="Fitted")
    }
    g <- g+ geom_line(aes(x = Time,y=Intensity,linetype=Model),data = s,size=1)
  }
  return(g)
}


plotLmms <- function(model,smooth,name,mean,...){
  Time <- Intensity <- Model <- NULL
  if(missing(mean))
    mean <- F
  if(missing(smooth))
    smooth <- F
  cl <- class(model)
  p <- NULL
  
  if(cl=="lm"){

    dfmain <- data.frame(Intensity=model$model$Expr,Time=model$model$time,size=1,Model="Mean")
    g <- ggplot() + geom_point(aes(x=Time,y=Intensity),data = dfmain,na.rm = T) +ggtitle(name)
    if(mean)
      g <- g + stat_summary(aes(x=Time,y=Intensity,linetype=Model),size=1,data=dfmain,fun.y=function(x)mean(x), geom="line",na.rm = T)
     
    if(smooth){ 
      spl <- spline(x = model$model$time, y = fitted(model), n = 500, method = "natural")
      s<- data.frame(Time=spl$x,Intensity=spl$y,Model="Smooth")
    }else{
      s <- data.frame(Intensity=fitted(model),Time=model$model$time,Model="Fitted")
    }
    
    g <- g+ geom_line(aes(x = Time,y=Intensity,linetype=Model),data = s,size=1)
    
  }else if(cl=='lme'){
    dfmain <- data.frame(Intensity=model$data$Expr,Time=model$data$time,size=1,Model="Mean")
    g <- ggplot()+ geom_point(aes(x=Time,y=Intensity),alpha=0.5,size=3,data = dfmain,na.rm = T) +ggtitle(name)
    if(mean)
      g <- g + stat_summary(aes(x=Time,y=Intensity,linetype=Model),size=1,data=dfmain,fun.y=function(x)mean(x), geom="line",na.rm = T)

   
    f <- fitted(model,level=1)
    if(smooth){
      spl <- spline(x = model$data$time, y = f, n = 500, method = "natural")
      s <- data.frame(Time=spl$x,Intensity=spl$y,Model="Smooth")
    }else{
      s <- data.frame(Intensity=na.omit(f),Time=model$data$time[!is.na(f)],Model="Fitted")
    }
  }else{
    return('Error plotting molecule')
  }
  g<-  g+ geom_line(aes(x = Time,y=Intensity,linetype=Model),size=1,data =s)
  return(g)
}

plotLmmsDeriv <- function(model,smooth,data,name,...){
  Time <- Intensity <- Model <- NULL
  if(missing(smooth))
    smooth <- F
  cl <- class(model)

  if(cl=="lm"){
    dfmain <- data.frame(Intensity=model$model$Expr,Time=model$model$time,Model="Mean")
    g <- ggplot()+ geom_point(aes(x=Time,y=Intensity),data = dfmain,na.rm = T) + ggtitle(name)
    g <- g + stat_summary(aes(x=Time,y=Intensity,linetype=Model),size=1,data=dfmain,fun.y=function(x)mean(x), geom="line",na.rm = T)

    if(smooth){
      spl <- spline(x =  as.numeric(colnames(data)), y = as.numeric(as.character(data)), n = 500, method = "natural")
      s <- data.frame(Time=spl$x,Intensity=spl$y,Model="Smooth")
    }else{
      s <- data.frame(Time=as.numeric(colnames(data)),Intensity=as.numeric(as.character(data)),Model="Derivative")
    }
    g <- g+geom_line(aes(x = Time,y=Intensity,linetype=Model),size=1,data =s)
  }else if(cl=='lme'){
    dfmain <- data.frame(Intensity=model$data$Expr,Time=model$data$time,size=1,Model="Mean")
    g <- ggplot()+ geom_point(aes(x=Time,y=Intensity),alpha=0.5,size=3,data = dfmain,na.rm = T) + ggtitle(name)
    
    g <- g + stat_summary(aes(x=Time,y=Intensity,linetype=Model),size=1,data=dfmain,fun.y=function(x)mean(x), geom="line",na.rm = T)

    if(smooth){
      spl <- spline(x =as.numeric(colnames(data)), y = data, n = 500, method = "natural")
      s <- data.frame(Time=spl$x,Intensity=spl$y,Model="Smooth")
    }else{
      s <- data.frame(Intensity=as.numeric(as.character(data)),Time=as.numeric(colnames(data)),Model="Derivative")
    }
    g<-  g+ geom_line(aes(x = Time,y=Intensity,linetype=Model),size=1,data =s)
  }else{
    g <- 'Error plotting molecule'
  }
  
  return(g)
  
}
