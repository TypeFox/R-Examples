
#Copyright (c) 2009-2014 Sebastien Bihorel
#All rights reserved.
#
#This file is part of scaRabee.
#
#    scaRabee is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    scaRabee is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with scaRabee.  If not, see <http://www.gnu.org/licenses/>.
#

estimation.plot <- function(problem=NULL,Fit=NULL,files=NULL){
  
# Set some default variable value
  nestparam <- length(Fit$estimations)
  
# Create ... vs iteration plots

  # Determine the most esthetical arrangement of plots based on nestparam
  tmplayout <- get.layout(nplot=nestparam + 1)
  
  # Import and reshape data to be plotted
  data <- read.table(file=files$iter,
                     header=TRUE,
                     as.is=TRUE,
                     sep=',',
                     check.names=FALSE)
  est.names <- problem$init$names[which(problem$init$isfix==0)]
  data <- reshape(data[,c('Objective function',est.names,'ID')],
                  idvar='iteration',
                  times=c('Objective function',est.names),
                  timevar='variable',
                  varying=list(c('Objective function',est.names)),
                  direction='long')
  data$variable <- factor(data$variable,levels=unique(data$variable))
  firstiterations <- data$iteration[match(unique(data$ID),data$ID)]
  for (id in unique(data$ID)){
    data$iteration[which(data$ID==id)] <- data$iteration[which(data$ID==id)] -
                                          firstiterations[id]
  }
  names(data)[3] <- 'value'
  
  # Objective function and parameters vs iteration plot
  pdf(file='01.Iterations.pdf',
      width=6.5,
      height=9)

  for (i in unique(data$ID)) {
    if (problem$method=='subject'){
      ititle <- paste('Subject', i)
    } else {
      ititle <- 'Population'
    }
    
    itplot <- xyplot(value~iteration|variable,
                     data=data,
                     subset=which(data$ID==i),
                     as.table=TRUE,
                     type='l',
                     col=4,
                     scales=list(x=list(relation='free'),
                                 y=list(relation='free')),
                     strip=strip.custom(style=1,bg=0),
                     xlab='Iteration',
                     ylab='Value / Estimate',
                     main=paste('Gradient monitoring - ',ititle,sep=''),
                     layout=tmplayout)
    print(itplot)
  }
  dev.off()
  
# Create diagnostic plots
  
  # Import data to be plotted
  data <- read.table(file=files$pred,
                     header=TRUE,
                     as.is=TRUE,
                     sep=',',
                     check.names=FALSE)
  
  # Creates prediction plots
  pdf(file='02.Predictions.pdf',
      width=6.5,
      height=9)
    
  for (iid in unique(data$ID)){
    subdata <- data[which(data$ID==iid),]
    for (itrt in unique(subdata$TRT)){
      if (problem$method=='subject'){
        ititle <- paste('Subject', iid, '- Treatment', itrt)
      } else {
        ititle <- paste('Population - Treatment', itrt)
      }
      trtdata <- subdata[which(subdata$TRT==itrt),]
      trtdata$DVID <- factor(trtdata$DVID,levels=unique(trtdata$DVID))
      noutput   <- nlevels(trtdata$DVID)
      tmplayout <- get.layout(nplot=noutput)
      
      myplot <- xyplot(DV+IPRED~TIME|DVID,
                       data=trtdata,
                       as.table=TRUE,
                       type=c('p','l'),
                       pch=c(3,NULL),
                       col=c(1,4),
                       distribute.type=TRUE,
                       scales=list(x=list(relation='free'),
                                   y=list(relation='free')),
                       strip=strip.custom(var.name='Output',
                                          style=1,
                                          bg=0,
                                          sep=' ',
                                          strip.name=c(TRUE,TRUE)),
                       xlab='Time',
                       ylab='Observations / Predictions',
                       main=ititle,
                       layout=tmplayout)
      print(myplot)
    }
  }
  
  dev.off()
  
  # Creates residuals plots
  pdf(file='03.Residuals.pdf',
      width=6.5,
      height=9)

  for (iid in unique(data$ID)){
    subdata <- data[which(data$ID==iid),]
    for (itrt in unique(subdata$TRT)){
      trtdata <- subdata[which(subdata$TRT==itrt),]
      for (icmt in unique(trtdata$DVID)){
        if (problem$method=='subject'){
          ititle <- paste('Subject', iid, '- Treatment', itrt, '- Output', icmt)
        } else {
          ititle <- paste('Population - Treatment', itrt, '- Output', icmt)
        }
        cmtdata <- trtdata[which(trtdata$DVID==icmt),]
        rplot1 <- xyplot(IPRED~DV,
                         data=cmtdata,
                         panel=function(x,y,...){
                           panel.xyplot(x,y,type='p',pch=3,col=4)
                           panel.abline(c(0,1),type='l',col=2)
                         },
                         strip=NULL,
                         xlab='Observations',
                         ylab='Predictions')
        rplot2 <- xyplot(WRES~TIME,
                         data=cmtdata,
                         panel=function(x,y,...){
                           panel.xyplot(x,y,type='p',pch=3,col=4)
                           panel.abline(c(0,0),type='l',col=2)
                         },
                         strip=NULL,
                         xlab='Time',
                         ylab='Weighted residuals')
        rplot3 <- xyplot(WRES~DV,
                         data=cmtdata,
                         panel=function(x,y,...){
                           panel.xyplot(x,y,type='p',pch=3,col=4)
                           panel.abline(c(0,0),type='l',col=2)
                         },
                         strip=NULL,
                         xlab='Observations',
                         ylab='Weighted residuals')
        rplot4 <- xyplot(WRES~IPRED,
                         data=cmtdata,
                         panel=function(x,y,...){
                           panel.xyplot(x,y,type='p',pch=3,col=4)
                           panel.abline(c(0,0),type='l',col=2)
                         },
                         strip=NULL,
                         xlab='Predictions',
                         ylab='Weighted residuals')
      
      print(rplot1,position=c(0,0.1,0.5,0.5),more=TRUE)
      grid.text(label=ititle,
                x=0.5,y=0.95,just='center',gp=gpar(cex=1.5))
      print(rplot2,position=c(0.5,0.1,1,0.5),more=TRUE)
      print(rplot3,position=c(0,0.5,0.5,0.9),more=TRUE)
      print(rplot4,position=c(0.5,0.5,1,0.9),more=FALSE)
      
      }
    }
  }
  dev.off()
  
  # Print a message to screen
  cat(sprintf('\n%s%s\n%s\n','Diagnostic figures have been created and saved in: ',
              getwd(),'You may open and edit them at your convenience.'))
                            
}
