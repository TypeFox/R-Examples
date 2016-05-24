# Jasmin Straube, Queensland Facility of Advanced Bioinformatics
# Part of this script was borrowed from the graphics and stats package.
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

#' Plot of \code{lmmsde} objects
#' 
#' Plot of the raw data the mean and the fitted \code{lmmsde} profile.
#' 
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom  stats spline na.omit
#' @param x An object of class \code{lmmsde}.
#' @param y \code{numeric} or \code{character} value. Either the row index or the row name determining which feature should be plotted. 
#' @param data alternative \code{matrix} or \code{data.frame} containing the original data for visualisation purposes.
#' @param time alternative \code{numeric} indicating the sample time point. Vector of same length as row lenghth of data for visualisation purposes.
#' @param group alternative \code{numeric} indicating the sample group. Vector of same length as row length of data for visualisation purposes.
#' @param type a \code{character} indicating what model to plot. Default  \code{'all'}, options: \code{'time'}, \code{'group'},\code{'group*time'}.
#' @param smooth an optional \code{logical} value. By default set to \code{FALSE}. If \code{TRUE} smooth representation of the fitted values. 
#' @param mean alternative \code{logical} if the mean should be displayed.  By default set to \code{TRUE}.
#' @param \ldots Additional arguments which are passed to \code{plot}.
#' @return plot showing raw data, mean profile and fitted profile. 
#' @examples 
#' \dontrun{
#' data(kidneySimTimeGroup)
#' lmmsDEtestl1 <-lmmsDE(data=kidneySimTimeGroup$data,time=kidneySimTimeGroup$time,
#'                 sampleID=kidneySimTimeGroup$sampleID,
#'                 group=kidneySimTimeGroup$group,
#'                 experiment="longitudinal1",basis="p-spline",keepModels=T) 
#' plot(lmmsDEtestl1,y=2,type="all")
#' plot(lmmsDEtestl1,y=2,type="time")
#' plot(lmmsDEtestl1,y=2,type="group")
#' plot(lmmsDEtestl1,y=2,type="group*time",smooth=TRUE)
#' 
#' #to save memory do not keep the models
#' lmmsDEtestl1 <-lmmsDE(data=kidneySimTimeGroup$data,time=kidneySimTimeGroup$time,
#'                 sampleID=kidneySimTimeGroup$sampleID,
#'                 group=kidneySimTimeGroup$group,
#'                 experiment="longitudinal1",basis="p-spline",keepModels=F) 
#'# just the fitted trajectory                 
#' plot(lmmsDEtestl1,y=2,type="all")
#' 
#' plot(lmmsDEtestl1,y=2,type="all",data=kidneySimTimeGroup$data,time=kidneySimTimeGroup$time,
#' group=kidneySimTimeGroup$group)}  

#' @method plot lmmsde
#' @export
plot.lmmsde <- function(x, y, data, time,group,type, smooth, mean,...){
 # library(graphics)

  if(missing(type)|sum(type%in%"all")>0){
    type <- c()
    if(length(x@modelGroup)>0|ncol(x@predGroup)>1)
      type <- c(type,'group')
    if(length(x@modelTime)>0|ncol(x@predTime)>1)
      type <- c(type,'time')
    if(length(x@modelTimeGroup)>0|ncol(x@predTimeGroup)>1)
      type <- c(type,'group*time')
  }
   name <- y
  if(length(grep("all",type))>0){
    type <- c('time','group','group*time') 
  }
  
  if(class(y)=='numeric')
    name <- as.character(x@DE$Molecule[y])
    
  if(class(y)=='character'|class(y)=="factor"){
    nam <- as.character(x@DE$Molecule)
    if(sum(nam%in%as.character(y))>0){
      name <- as.character(y)
      y <-which(nam%in%y)
    }else{
      stop(paste('Could not find feature',y,'in rownames(x@pred.spline).'))
    }
  }

  if(missing(mean))
    mean <- F
  
  if(sum(type%in%'time')>0){
    name2 <- paste(name,'time')
    if(length(x@modelTime)>0){
     p1 <- plotLmms(x@modelTime[[y]],smooth=smooth,name2,mean=mean,...)

    }else{
      p1 <- plotModel(x@predTime,y,smooth=smooth,name2,data=data,time2=time,mean=mean)
    }
    if(length(type)<3){
      suppressWarnings(print(p1))
    }
  }
  if(sum(type%in%"group")>0){
   name2 <- paste(name,'group')
   if(length(x@modelGroup)>0){
      p2 <- plotLmmsdeFunc(x@modelGroup,index=y,smooth=smooth,name2,mean,...)
   }else{
     p2 <- plotModel(x@predGroup,y,smooth=smooth,name2,data=data,time2=time,group=group,mean=mean)
   }
   
   if(length(type)<3){
     suppressWarnings(print(p2))
   }
   
  }
  if(sum(type%in%"group*time")>0){
    name2 <- paste(name,'group*time')
    if(length(x@modelTimeGroup)>0){
      p3 <- plotLmmsdeFunc(x@modelTimeGroup,index=y,smooth=smooth,name2,mean,...)
    } else{
     p3 <- plotModel(x@predTimeGroup,index = y,smooth=smooth,name2,data=data,time2=time,group=group,mean=mean)
    }
    if(length(type)<3){
      suppressWarnings(print(p3))
    }
  }
  
  if(length(type)==3|sum(type%in%"all")>0){
    grid.arrange(p1, p2, p3, ncol=2)
  }
  }


plotLmmsdeFunc <- function(object,index,smooth,name,mean,...){
  Time <- Intensity <- Model <-Group<-  NULL
  if(missing(smooth))
    smooth <- F
  model <- object[[index]]

  if(is.null(model))
    stop("Requested model not available")
  cl <- class(model)
    
    group <- model$data$Group
    g1 <-which(group==unique(group)[1])
   
    g2 <- which(group==unique(group)[2])

    g1Label <-sort(unique(group))[1]
    g2Label <- sort(unique(group))[2]
  
    dfmain <- data.frame(Intensity=model$data$Expr,Time=model$data$time,Group=group,size=1)
    g <- ggplot()+ geom_point(aes(x=Time,y=Intensity,shape=Group,color=Group),alpha=0.5,size=3,data = dfmain,na.rm = T) +ggtitle(name)
    if(mean)
      g <- g + stat_summary(aes(x=Time,y=Intensity,group=Group,colour=Group,linetype='Mean'),size=1,data=dfmain,fun.y=function(x)mean(x), geom="line",na.rm = T)

    f <- fitted(model,level=1)
    f1 <- f[g1]
    f2<- f[g2]
  
    dfmain <- data.frame(Intensity=f,Time=model$data$time,Group=group,size=1,Model="Fitted")
  
  
  
    if(smooth){
      spl1 <- spline(x = model$data$time[g1], y = f1, n = 500, method = "natural")
      s1 <- data.frame(Time=spl1$x,Intensity=spl1$y,Model="Smooth",Group=g1Label)
      
      spl2 <- spline(x = model$data$time[g2], y = f2, n = 500, method = "natural")
      s2 <- data.frame(Time=spl2$x,Intensity=spl2$y,Model="Smooth",Group=g2Label)
    }else{
      s1 <- data.frame(Intensity=na.omit(f1),Time=model$data$time[intersect(which(!is.na(f)),g1)],Model="Fitted",Group=g1Label)
      
      s2 <- data.frame(Intensity=na.omit(f2),Time=model$data$time[intersect(which(!is.na(f)),g2)],Model="Fitted",Group=g2Label)
    }
    g <- g+ geom_line(aes(x = Time,y=Intensity,color=Group,linetype=Model),data = s1[!is.na(s1$Intensity),],size=1)+ geom_line(aes(x = Time,y=Intensity,color=Group,linetype=Model),size=1,data = s2[!is.na(s2$Intensity),])
  
  return(g)
 }

plotModel <- function(object,index,smooth,name,data,time2,group,mean,...){  
  Time <- Intensity <- Model <- Group<- NULL
  if(missing(smooth))
    smooth <- F
  if(missing(mean))
    mean <- F

  model <- object[index,]
 if(sum(is.na(model))==length(model))
    stop('Error plotting molecule')
  time <- suppressWarnings(as.numeric(names(model)))

  s <- strsplit(colnames(object),split = " ")
  s <- sapply(s,'[')

  
  if(!is.null(dim(s)) & !is.null(s)){
    group2 <- s[1,]
   
    time <-  suppressWarnings(as.numeric(s[2,]))
    g1Label <-na.omit(sort(unique(group2)))[1]
    g2Label <- na.omit(sort(unique(group2)))[2]
    
    g1 <- which(group2==g1Label)
    g2 <- which(group2==g2Label)
    if(!missing(group)){
      group2[g1] <- unique(group)[1]
      group2[g2] <- unique(group)[2]
    }
    if(!missing(data) & !missing(time2) & !missing(group)){
      dfmain <<- data.frame(Intensity=data[,index],Time=time2,Group=factor(group),size=1,Model="Mean")
      g <- ggplot()+ geom_point(aes(x=Time,y=Intensity,shape=Group,color=Group),alpha=0.5,size=3,data = dfmain,na.rm = T) + ggtitle(name)
      if(mean)
        g <- g + stat_summary(aes(x=Time,y=Intensity,linetype=Model,color=Group,group=Group),size=1,data=dfmain,fun.y=function(x)mean(x), geom="line",na.rm = T)
    }else{
      dfmain <- data.frame(Intensity=model,Time=time,Group=factor(group2),size=1,Model="Mean")
      g <- ggplot()+ geom_point(aes(x=Time,y=Intensity,shape=Group,color=Group),alpha=0.5,size=3,data = dfmain,na.rm = T) + ggtitle(name)
      if(mean)
        g <- g + stat_summary(aes(x=Time,y=Intensity,linetype=Model,color=Group,group=Group),size=1,data=dfmain,fun.y=function(x)mean(x), geom="line",na.rm = T)
    
    }
    if(smooth){
      spl1 <- as.data.frame(spline(x = time[g1], y = model[g1], n = 500, method = "natural"))
      s1 <- data.frame(Time=spl1$x,Intensity=spl1$y,Model="Smooth",Group=g1Label)
      
      spl2 <- as.data.frame(spline(x = time[g2], y = model[g2], n = 500, method = "natural"))
      s2 <- data.frame(Time=spl2$x,Intensity=spl2$y,Model="Smooth",Group=g2Label)
     
    }else{
      s1 <- data.frame(Time = time[g1], Intensity = model[g1],Model="Fitted",Group=g1Label)
      s2 <- data.frame(Time = time[g2], Intensity = model[g2],Model="Fitted",Group=g2Label)
    }
    g <- g+ geom_line(aes(x = Time,y=Intensity,color=Group,linetype=Model),data = s1[!is.na(s1$Intensity),])+ geom_line(aes(x = Time,y=Intensity,color=Group,linetype=Model),data = s2[!is.na(s2$Intensity),],size=1)    
     
    }else{
   # time <- as.numeric(as.character(s))
    if(!missing(data) & !missing(time2)){
      dfmain <- data.frame(Intensity=data[,index],Time=time2)
    }else{
      dfmain <- data.frame(Intensity=model,Time=time)    
    }
    
    g <- ggplot()+ geom_point(aes(x=Time,y=Intensity),alpha=0.5,size=3,data = dfmain,na.rm = T) +ggtitle(name)
    
    if(smooth){
      sp1 <- as.data.frame(spline(x = time, y = model, n = 500, method = "natural"))
      dfModel <- data.frame(Intensity=sp1$y,Time=sp1$x,Model="Smooth")                        
    }else{
      dfModel <- data.frame(Intensity=model,Time=time,Model="Fitted")
    }

     g<-g+ geom_line(aes(x = Time,y=Intensity,linetype=Model),data =dfModel[!is.na(dfModel$Intensity),],size=1)
     if(mean)
       g <- g + stat_summary(aes(x=Time,y=Intensity,linetype='Mean'),size=1,data=dfmain,fun.y=function(x)mean(x), geom="line",na.rm = T)
     
     }
                     
  return(g)
  

}
