plot.GroAgeData <-
function(x, ...)
    {
     options(warn=-1)
     spec.name <- x$Properties$Species
     #One group
     if(x$Properties$Sex == "Females" | x$Properties$Sex == "Males" | x$Properties$Sex == "Pooled")
       {
        if(is.na(x$Properties["Units.Weight"]) & is.na(x$Properties["Units.Liver"]) & is.na(x$Properties["Units.Gonad"]))
          {
           age    <- x$Data[,2]
           len    <- x$Data[,3]
           par(mfrow=c(3,1),oma=c(2,2,2,1), mar=c(4,4,2,2))
           hist(x=age,main="",xlab=paste("Age (",x$Properties["Units.Age"],")",sep=""),ylab="Frequency")
           hist(x=len,main="",xlab=paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),ylab="Frequency")
           plot(x=age,y=len,xlab=paste("Age (",x$Properties["Units.Age"],")",sep=""),ylab=paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),pch=19)
           mtext(side=3,text=x$Properties$Sex,outer=TRUE)
          }
        else if(!is.na(x$Properties["Units.Weight"]) & is.na(x$Properties["Units.Gonad"]) & is.na(x$Properties["Units.Liver"]))
          {
           age    <- x$Data[,2]
           len    <- x$Data[,3]
           wgt    <- x$Data[,4]
           par(mfrow=c(2,2),oma=c(2,2,2,1), mar=c(4,4,2,2))
           hist(x=age,main="",xlab=paste("Age (",x$Properties["Units.Age"],")",sep=""),ylab="Frequency")
           hist(x=len,main="",xlab=paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),ylab="Frequency")
           plot(x=age,y=len,xlab=paste("Age (",x$Properties["Units.Age"],")",sep=""),ylab=paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),pch=19)
           plot(x=len,y=wgt,xlab=paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),ylab=paste("Weight (",x$Properties["Units.Weight"],")",sep=""),pch=19)
           mtext(side=3,text=x$Properties$Sex,outer=TRUE)
          }
        else if(!is.na(x$Properties["Units.Weight"]) & !is.na(x$Properties["Units.Gonad"]) & !is.na(x$Properties["Start.Month"]) & is.na(x$Properties["Units.Liver"]))
          {
           age    <- x$Data[,2]
           len    <- x$Data[,3]
           wgt    <- x$Data[,4]
           if(x$Properties["Start.Month"] == x$Properties["End.Month"])
             {
              par(mfrow=c(2,2),oma=c(2,2,2,1), mar=c(4,4,2,2))
              hist(x=age,main="",xlab=paste("Age (",x$Properties["Units.Age"],")",sep=""),ylab="Frequency")
              hist(x=len,main="",xlab=paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),ylab="Frequency")
              plot(x=age,y=len,xlab=paste("Age (",x$Properties["Units.Age"],")",sep=""),ylab=paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),pch=19)
              plot(x=len,y=wgt,xlab=paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),ylab=paste("Weight (",x$Properties["Units.Weight"],")",sep=""),pch=19)
              mtext(side=3,text=x$Properties$Sex,outer=TRUE)
             }
           else
             {
              gsi <- 100*x$Data[,5]/wgt
              mon <- as.numeric(format(x$Data$Date, "%m"))
              layout(matrix(c(1,2,3,4,5,5),3,2, byrow = TRUE))
              par(oma=c(2,2,2,1), mar=c(4,4,2,2))
              hist(x=age,main="",xlab=paste("Age (",x$Properties["Units.Age"],")",sep=""),ylab="Frequency")
              hist(x=len,main="",xlab=paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),ylab="Frequency")
              plot(x=age,y=len,xlab=paste("Age (",x$Properties["Units.Age"],")",sep=""),ylab=paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),pch=19)
              plot(x=len,y=wgt,xlab=paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),ylab=paste("Weight (",x$Properties["Units.Weight"],")",sep=""),pch=19)
              plot(1:12,rep(0,12),ylim=c(0,max(gsi,na.rm=TRUE)),ylab="GSI (%)", xlab="Month",type="n")
              points(x=mon,y=gsi, pch=19)
              mtext(side=3,text=x$Properties$Sex,outer=TRUE)
             }
          }
        else if(!is.na(x$Properties["Units.Weight"]) & !is.na(x$Properties["Units.Gonad"]) & !is.na(x$Properties["Start.Month"]) & !is.na(x$Properties["Units.Liver"]))
          {
           age    <- x$Data[,2]
           len    <- x$Data[,3]
           wgt    <- x$Data[,4]
           if(x$Properties["Start.Month"] == x$Properties["End.Month"])
             {
              par(mfrow=c(2,2),oma=c(2,2,2,1), mar=c(4,4,2,2))
              hist(x=age,main="",xlab=paste("Age (",x$Properties["Units.Age"],")",sep=""),ylab="Frequency")
              hist(x=len,main="",xlab=paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),ylab="Frequency")
              plot(x=age,y=len,xlab=paste("Age (",x$Properties["Units.Age"],")",sep=""),ylab=paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),pch=19)
              plot(x=len,y=wgt,xlab=paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),ylab=paste("Weight (",x$Properties["Units.Weight"],")",sep=""),pch=19)
              mtext(side=3,text=x$Properties$Sex,outer=TRUE)
             }
           else
             {
              gsi <- 100*x$Data[,5]/wgt
              hsi <- 100*x$Data[,6]/wgt
              mon <- as.numeric(format(x$Data$Date, "%m"))
              layout(matrix(c(1,2,3,4,5,5,6,6),4,2, byrow = TRUE))
              par(oma=c(2,2,2,1), mar=c(4,4,2,2))
              hist(x=age,main="",xlab=paste("Age (",x$Properties["Units.Age"],")",sep=""),ylab="Frequency")
              hist(x=len,main="",xlab=paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),ylab="Frequency")
              plot(x=age,y=len,xlab=paste("Age (",x$Properties["Units.Age"],")",sep=""),ylab=paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),pch=19)
              plot(x=len,y=wgt,xlab=paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),ylab=paste("Weight (",x$Properties["Units.Weight"],")",sep=""),pch=19)
              plot(1:12,rep(0,12),ylim=c(0,max(gsi,na.rm=TRUE)),ylab="GSI (%)", xlab="Month",type="n")
              points(x=mon,y=gsi, pch=19)
              plot(1:12,rep(0,12),ylim=c(0,max(hsi,na.rm=TRUE)),ylab="HSI (%)", xlab="Month",type="n")
              points(x=mon,y=hsi, pch=19)
              mtext(side=3,text=x$Properties$Sex,outer=TRUE)
             }
          }
       }
     #Two groups (males, females)
     else if(x$Properties$Sex == "Both")
       {
        if(is.na(x$Properties["Units.Weight"]) & is.na(x$Properties["Units.Liver"]) & is.na(x$Properties["Units.Gonad"]))
          {
           age     <- c(x$Data$Females[,2],x$Data$Males[,2])
           age.tag <- c(rep("Females",length(x$Data$Females[,2])),rep("Males",length(x$Data$Males[,2])))
           len     <- c(x$Data$Females[,3],x$Data$Males[,3])
           len.tag <- c(rep("Females",length(x$Data$Females[,3])),rep("Males",length(x$Data$Males[,3])))
           par(mfrow=c(3,1),oma=c(2,2,2,1), mar=c(4,4,2,2))
           Ecdf(age,group=age.tag,xlab="",
                ylab="",col=c("red","blue"),label.curves=FALSE)
           mtext(expression("Proportion " <= " age"),2,line=2)
           mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2)
           Ecdf(len,group=len.tag,xlab="",
                ylab="",col=c("red","blue"),label.curves=FALSE)
           mtext(expression("Proportion " <= " length"),2,line=2)
           mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2)
           plot(x=age[1:length(x$Data$Females[,2])],y=len[1:length(x$Data$Females[,3])],
                xlab="",ylab="",pch=19,
                xlim=c(min(age,na.rm=TRUE),max(age,na.rm=TRUE)),ylim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),col="red")
           mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2)
           mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),2,line=2)
           points(x=age[(length(x$Data$Females[,2])+1):(length(x$Data$Females[,2])+length(x$Data$Males[,2]))],
                  y=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                  pch=19,col="blue")
           mtext(side=3,text="Females",outer=TRUE,col="red",line=-1)
           mtext(side=3,text="Males",outer=TRUE,col="blue",line=0)
          }
        else if(!is.null(x$Properties["Units.Weight"]) & is.null(x$Properties["Units.Gonad"]) & is.null(x$Properties["Units.Liver"]))
          {
           age     <- c(x$Data$Females[,2],x$Data$Males[,2])
           age.tag <- c(rep("Females",length(x$Data$Females[,2])),rep("Males",length(x$Data$Males[,2])))
           len     <- c(x$Data$Females[,3],x$Data$Males[,3])
           len.tag <- c(rep("Females",length(x$Data$Females[,3])),rep("Males",length(x$Data$Males[,3])))
           wgt     <- c(x$Data$Females[,4],x$Data$Males[,4])
           wgt.tag <- c(rep("Females",length(x$Data$Females[,4])),rep("Males",length(x$Data$Males[,4])))
           par(mfrow=c(2,2),oma=c(2,2,2,1), mar=c(4,4,2,2))
           Ecdf(age,group=age.tag,xlab="",
                ylab="",col=c("red","blue"),label.curves=FALSE)
           mtext(expression("Proportion " <= " age"),2,line=2)
           mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2)
           Ecdf(len,group=len.tag,xlab="",
                ylab="",col=c("red","blue"),label.curves=FALSE)
           mtext(expression("Proportion " <= " length"),2,line=2)
           mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2)
           plot(x=len[1:length(x$Data$Females[,3])],y=wgt[1:length(x$Data$Females[,4])],
                xlab="",ylab="",pch=19,
                xlim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),ylim=c(min(wgt,na.rm=TRUE),max(wgt,na.rm=TRUE)),col="red")
           mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2)
           mtext(paste("Body Weight (",x$Properties["Units.Weight"],")",sep=""),2,line=2)
           points(x=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                  y=wgt[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                  pch=19,col="blue")
           plot(x=age[1:length(x$Data$Females[,2])],y=len[1:length(x$Data$Females[,3])],
                xlab="",ylab="",pch=19,
                xlim=c(min(age,na.rm=TRUE),max(age,na.rm=TRUE)),ylim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),col="red")
           mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2)
           mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),2,line=2)
           points(x=age[(length(x$Data$Females[,2])+1):(length(x$Data$Females[,2])+length(x$Data$Males[,2]))],
                  y=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                  pch=19,col="blue")
           mtext(side=3,text="Females",outer=TRUE,col="red",line=-1)
           mtext(side=3,text="Males",outer=TRUE,col="blue",line=0)
          }
        else if(!is.null(x$Properties["Units.Weight"]) & !is.null(x$Properties["Units.Gonad"]) & !is.null(x$Properties["Start.Month"]) & is.null(x$Properties["Units.Liver"]))
          {
           age     <- c(x$Data$Females[,2],x$Data$Males[,2])
           age.tag <- c(rep("Females",length(x$Data$Females[,2])),rep("Males",length(x$Data$Males[,2])))
           len     <- c(x$Data$Females[,3],x$Data$Males[,3])
           len.tag <- c(rep("Females",length(x$Data$Females[,3])),rep("Males",length(x$Data$Males[,3])))
           wgt     <- c(x$Data$Females[,4],x$Data$Males[,4])
           wgt.tag <- c(rep("Females",length(x$Data$Females[,4])),rep("Males",length(x$Data$Males[,4])))
           if(x$Properties["Start.Month"] == x$Properties["End.Month"])
             {
              par(mfrow=c(2,2),oma=c(2,2,2,1), mar=c(4,4,2,2))
              Ecdf(age,group=age.tag,xlab="",
                   ylab="",col=c("red","blue"),label.curves=FALSE)
              mtext(expression("Proportion " <= " age"),2,line=2)
              mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2)
              Ecdf(len,group=len.tag,xlab="",
                   ylab="",col=c("red","blue"),label.curves=FALSE)
              mtext(expression("Proportion " <= " length"),2,line=2)
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2)
              plot(x=len[1:length(x$Data$Females[,3])],y=wgt[1:length(x$Data$Females[,4])],
                   xlab="",ylab="",pch=19,
                   xlim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),ylim=c(min(wgt,na.rm=TRUE),max(wgt,na.rm=TRUE)),col="red")
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2)
              mtext(paste("Body Weight (",x$Properties["Units.Weight"],")",sep=""),2,line=2)
              points(x=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                     y=wgt[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                     pch=19,col="blue")
              plot(x=age[1:length(x$Data$Females[,2])],y=len[1:length(x$Data$Females[,3])],
                   xlab="",ylab="",pch=19,
                   xlim=c(min(age,na.rm=TRUE),max(age,na.rm=TRUE)),ylim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),col="red")
              mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2)
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),2,line=2)
              points(x=age[(length(x$Data$Females[,2])+1):(length(x$Data$Females[,2])+length(x$Data$Males[,2]))],
                     y=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                     pch=19,col="blue")
              mtext(side=3,text="Females",outer=TRUE,col="red",line=-1)
              mtext(side=3,text="Males",outer=TRUE,col="blue",line=0)
             }
           else
             {
              gsi     <- c(100*x$Data$Females[,5]/x$Data$Females[,4],100*x$Data$Males[,5]/x$Data$Males[,4])
              gsi.tag <- c(rep("Females",length(x$Data$Females[,5])),rep("Males",length(x$Data$Males[,5])))
              mon     <- c(as.numeric(format(x$Data$Females$Date, "%m")),as.numeric(format(x$Data$Males$Date, "%m")))
              mon.tag <- c(rep("Females",length(x$Data$Females$Date)),rep("Males",length(x$Data$Males$Date)))
              layout(matrix(c(1,2,3,4,5,5),3,2, byrow = TRUE))
              par(oma=c(2,2,2,1), mar=c(4,4,2,2))
              Ecdf(age,group=age.tag,xlab="",
                   ylab="",col=c("red","blue"),label.curves=FALSE)
              mtext(expression("Proportion " <= " age"),2,line=2)
              mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2)
              Ecdf(len,group=len.tag,xlab="",
                   ylab="",col=c("red","blue"),label.curves=FALSE)
              mtext(expression("Proportion " <= " length"),2,line=2)
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2)
              plot(x=len[1:length(x$Data$Females[,3])],y=wgt[1:length(x$Data$Females[,4])],
                   xlab="",ylab="",pch=19,
                   xlim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),ylim=c(min(wgt,na.rm=TRUE),max(wgt,na.rm=TRUE)),col="red")
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2)
              mtext(paste("Body Weight (",x$Properties["Units.Weight"],")",sep=""),2,line=2)
              points(x=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                     y=wgt[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                     pch=19,col="blue")
              plot(x=age[1:length(x$Data$Females[,2])],y=len[1:length(x$Data$Females[,3])],
                   xlab="",ylab="",pch=19,
                   xlim=c(min(age,na.rm=TRUE),max(age,na.rm=TRUE)),ylim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),col="red")
              mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2)
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),2,line=2)
              points(x=age[(length(x$Data$Females[,2])+1):(length(x$Data$Females[,2])+length(x$Data$Males[,2]))],
                     y=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                     pch=19,col="blue")
              plot(1:12,rep(0,12),ylim=c(0,max(gsi,na.rm=TRUE)),ylab="", xlab="",type="n")
              mtext("GSI (%)",2,line=2)
              mtext("Month",1,line=2)
              points(x=mon[1:length(x$Data$Females[,4])],y=gsi[1:length(x$Data$Females[,4])], pch=19, col="red")
              points(x=mon[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                     y=gsi[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                     pch=19, col="blue")
              mtext(side=3,text="Females",outer=TRUE,col="red",line=-1)
              mtext(side=3,text="Males",outer=TRUE,col="blue",line=0)
             }
          }
        else if(!is.null(x$Properties["Units.Weight"]) & !is.null(x$Properties["Units.Gonad"]) & !is.null(x$Properties["Start.Month"]) & !is.null(x$Properties["Units.Liver"]))
          {
           age     <- c(x$Data$Females[,2],x$Data$Males[,2])
           age.tag <- c(rep("Females",length(x$Data$Females[,2])),rep("Males",length(x$Data$Males[,2])))
           len     <- c(x$Data$Females[,3],x$Data$Males[,3])
           len.tag <- c(rep("Females",length(x$Data$Females[,3])),rep("Males",length(x$Data$Males[,3])))
           wgt     <- c(x$Data$Females[,4],x$Data$Males[,4])
           wgt.tag <- c(rep("Females",length(x$Data$Females[,4])),rep("Males",length(x$Data$Males[,4])))
           if(x$Properties["Start.Month"] == x$Properties["End.Month"])
             {
              par(mfrow=c(2,2),oma=c(2,2,2,1), mar=c(4,4,2,2))
              Ecdf(age,group=age.tag,xlab="",
                   ylab="",col=c("red","blue"),label.curves=FALSE)
              mtext(expression("Proportion " <= " age"),2,line=2)
              mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2)
              Ecdf(len,group=len.tag,xlab="",
                   ylab="",col=c("red","blue"),label.curves=FALSE)
              mtext(expression("Proportion " <= " length"),2,line=2)
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2)
              plot(x=len[1:length(x$Data$Females[,3])],y=wgt[1:length(x$Data$Females[,4])],
                   xlab="",ylab="",pch=19,
                   xlim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),ylim=c(min(wgt,na.rm=TRUE),max(wgt,na.rm=TRUE)),col="red")
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2)
              mtext(paste("Body Weight (",x$Properties["Units.Weight"],")",sep=""),2,line=2)
              points(x=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                     y=wgt[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                     pch=19,col="blue")
              plot(x=age[1:length(x$Data$Females[,2])],y=len[1:length(x$Data$Females[,3])],
                   xlab="",ylab="",pch=19,
                   xlim=c(min(age,na.rm=TRUE),max(age,na.rm=TRUE)),ylim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),col="red")
              mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2)
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),2,line=2)
              points(x=age[(length(x$Data$Females[,2])+1):(length(x$Data$Females[,2])+length(x$Data$Males[,2]))],
                     y=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                     pch=19,col="blue")
              mtext(side=3,text="Females",outer=TRUE,col="red",line=-1)
              mtext(side=3,text="Males",outer=TRUE,col="blue",line=0)
             }
           else
             {
              gsi     <- c(100*x$Data$Females[,5]/x$Data$Females[,4],100*x$Data$Males[,5]/x$Data$Males[,4])
              gsi.tag <- c(rep("Females",length(x$Data$Females[,5])),rep("Males",length(x$Data$Males[,5])))
              hsi     <- c(100*x$Data$Females[,6]/x$Data$Females[,4],100*x$Data$Males[,6]/x$Data$Males[,4])
              hsi.tag <- c(rep("Females",length(x$Data$Females[,6])),rep("Males",length(x$Data$Males[,6])))
              mon     <- c(as.numeric(format(x$Data$Females$Date, "%m")),as.numeric(format(x$Data$Males$Date, "%m")))
              mon.tag <- c(rep("Females",length(x$Data$Females$Date)),rep("Males",length(x$Data$Males$Date)))
              layout(matrix(c(1,2,3,4,5,5,6,6),4,2, byrow = TRUE))
              par(oma=c(2,2,2,1), mar=c(4,4,2,2))
              Ecdf(age,group=age.tag,xlab="",
                   ylab="",col=c("red","blue"),label.curves=FALSE)
              mtext(expression("Proportion " <= " age"),2,line=2,cex=0.75)
              mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2,cex=0.75)
              Ecdf(len,group=len.tag,xlab="",
                   ylab="",col=c("red","blue"),label.curves=FALSE)
              mtext(expression("Proportion " <= " length"),2,line=2,cex=0.75)
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2,cex=0.75)
              plot(x=len[1:length(x$Data$Females[,3])],y=wgt[1:length(x$Data$Females[,4])],
                   xlab="",ylab="",pch=19,
                   xlim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),ylim=c(min(wgt,na.rm=TRUE),max(wgt,na.rm=TRUE)),col="red")
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2,cex=0.75)
              mtext(paste("Body Weight (",x$Properties["Units.Weight"],")",sep=""),2,line=2,cex=0.75)
              points(x=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                     y=wgt[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                     pch=19,col="blue")
              plot(x=age[1:length(x$Data$Females[,2])],y=len[1:length(x$Data$Females[,3])],
                   xlab="",ylab="",pch=19,
                   xlim=c(min(age,na.rm=TRUE),max(age,na.rm=TRUE)),ylim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),col="red")
              mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2,cex=0.75)
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),2,line=2,cex=0.75)
              points(x=age[(length(x$Data$Females[,2])+1):(length(x$Data$Females[,2])+length(x$Data$Males[,2]))],
                     y=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                     pch=19,col="blue")
              plot(1:12,rep(0,12),ylim=c(0,max(gsi,na.rm=TRUE)),ylab="", xlab="",type="n")
              mtext("GSI (%)",2,line=2,cex=0.75)
              mtext("Month",1,line=2,cex=0.75)
              points(x=mon[1:length(x$Data$Females[,4])],y=gsi[1:length(x$Data$Females[,4])], pch=19, col="red")
              points(x=mon[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                     y=gsi[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                     pch=19, col="blue")
              plot(1:12,rep(0,12),ylim=c(0,max(hsi,na.rm=TRUE)),ylab="", xlab="",type="n")
              mtext("HSI (%)",2,line=2,cex=0.75)
              mtext("Month",1,line=2,cex=0.75)
              points(x=mon[1:length(x$Data$Females[,4])],y=hsi[1:length(x$Data$Females[,4])], pch=19, col="red")
              points(x=mon[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                     y=hsi[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                     pch=19, col="blue")
              mtext(side=3,text="Females",outer=TRUE,col="red",line=-1)
              mtext(side=3,text="Males",outer=TRUE,col="blue",line=0)
             }
          }

       }
     #Three groups (Males, females, and unsexed)
     else
       {
        if(is.na(x$Properties["Units.Weight"]) & is.na(x$Properties["Units.Liver"]) & is.na(x$Properties["Units.Gonad"]))
          {
           age     <- c(x$Data$Females[,2],x$Data$Males[,2],x$Data$Unsexed[,2])
           age.tag <- c(rep("Females",length(x$Data$Females[,2])),rep("Males",length(x$Data$Males[,2])),rep("Unsexed",length(x$Data$Unsexed[,2])))
           len     <- c(x$Data$Females[,3],x$Data$Males[,3],x$Data$Unsexed[,3])
           len.tag <- c(rep("Females",length(x$Data$Females[,3])),rep("Males",length(x$Data$Males[,3])),rep("Unsexed",length(x$Data$Unsexed[,3])))
           par(mfrow=c(3,1),oma=c(2,2,2,1), mar=c(4,4,2,2))
           Ecdf(age,group=age.tag,xlab="",
                ylab="",col=c("red","blue","black"),label.curves=FALSE)
           mtext(expression("Proportion " <= " age"),2,line=2)
           mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2)
           Ecdf(len,group=len.tag,xlab="",
                ylab="",col=c("red","blue","black"),label.curves=FALSE)
           mtext(expression("Proportion " <= " length"),2,line=2)
           mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2)
           plot(x=age[1:length(x$Data$Females[,2])],y=len[1:length(x$Data$Females[,3])],
                xlab="",ylab="",pch=1,
                xlim=c(min(age,na.rm=TRUE),max(age,na.rm=TRUE)),ylim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),col="red")
           mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2)
           mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),2,line=2)
           points(x=age[(length(x$Data$Females[,2])+1):(length(x$Data$Females[,2])+length(x$Data$Males[,2]))],
                  y=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                  pch=1,col="blue")
           points(x=age[(length(x$Data$Females[,2])+length(x$Data$Males[,2])+1):(length(x$Data$Females[,2])+length(x$Data$Males[,2])+length(x$Data$Unsexed[,2]))],
                  y=len[(length(x$Data$Females[,3])+length(x$Data$Males[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3])+length(x$Data$Unsexed[,3]))],
                  pch=1,col="black")
           mtext(side=3,text="Unsexed",outer=TRUE,col="black",line=-2)
           mtext(side=3,text="Females",outer=TRUE,col="red",line=-1)
           mtext(side=3,text="Males",outer=TRUE,col="blue",line=0)
          }
        else if(!is.na(x$Properties["Units.Weight"]) & is.na(x$Properties["Units.Gonad"]) & is.na(x$Properties["Units.Liver"]))
          {
           age     <- c(x$Data$Females[,2],x$Data$Males[,2],x$Data$Unsexed[,2])
           age.tag <- c(rep("Females",length(x$Data$Females[,2])),rep("Males",length(x$Data$Males[,2])),rep("Unsexed",length(x$Data$Unsexed[,2])))
           len     <- c(x$Data$Females[,3],x$Data$Males[,3],x$Data$Unsexed[,3])
           len.tag <- c(rep("Females",length(x$Data$Females[,3])),rep("Males",length(x$Data$Males[,3])),rep("Unsexed",length(x$Data$Unsexed[,3])))
           wgt     <- c(x$Data$Females[,4],x$Data$Males[,4],x$Data$Unsexed[,4])
           wgt.tag <- c(rep("Females",length(x$Data$Females[,4])),rep("Males",length(x$Data$Males[,4])),rep("Unsexed",length(x$Data$Unsexed[,4])))
           par(mfrow=c(2,2),oma=c(2,2,2,1), mar=c(4,4,2,2))
           Ecdf(age,group=age.tag,xlab="",
                ylab="",col=c("red","blue","black"),label.curves=FALSE)
           mtext(expression("Proportion " <= " age"),2,line=2)
           mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2)
           Ecdf(len,group=len.tag,xlab="",
                ylab="",col=c("red","blue","black"),label.curves=FALSE)
           mtext(expression("Proportion " <= " length"),2,line=2)
           mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2)
           plot(x=len[1:length(x$Data$Females[,3])],y=wgt[1:length(x$Data$Females[,4])],
                xlab="",ylab="",pch=1,
                xlim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),ylim=c(min(wgt,na.rm=TRUE),max(wgt,na.rm=TRUE)),col="red")
           mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2)
           mtext(paste("Body Weight (",x$Properties["Units.Weight"],")",sep=""),2,line=2)
           points(x=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                  y=wgt[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                  pch=1,col="blue")
           points(x=len[(length(x$Data$Females[,3])+length(x$Data$Males[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3])+length(x$Data$Unsexed[,3]))],
                  y=wgt[(length(x$Data$Females[,4])+length(x$Data$Males[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4])+length(x$Data$Unsexed[,4]))],
                  pch=1,col="black")
           plot(x=age[1:length(x$Data$Females[,2])],y=len[1:length(x$Data$Females[,3])],
                xlab="",ylab="",pch=1,
                xlim=c(min(age,na.rm=TRUE),max(age,na.rm=TRUE)),ylim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),col="red")
           mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2)
           mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),2,line=2)
           points(x=age[(length(x$Data$Females[,2])+1):(length(x$Data$Females[,2])+length(x$Data$Males[,2]))],
                  y=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                  pch=1,col="blue")
           points(x=age[(length(x$Data$Females[,2])+length(x$Data$Males[,2])+1):(length(x$Data$Females[,2])+length(x$Data$Males[,2])+length(x$Data$Unsexed[,2]))],
                  y=len[(length(x$Data$Females[,3])+length(x$Data$Males[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3])+length(x$Data$Unsexed[,3]))],
                  pch=1,col="black")
           mtext(side=3,text="Unsexed",outer=TRUE,col="black",line=-2)
           mtext(side=3,text="Females",outer=TRUE,col="red",line=-1)
           mtext(side=3,text="Males",outer=TRUE,col="blue",line=0)
          }
        else if(!is.na(x$Properties["Units.Weight"]) & !is.na(x$Properties["Units.Gonad"]) & !is.na(x$Properties["Start.Month"]) & is.na(x$Properties["Units.Liver"]))
          {
           age     <- c(x$Data$Females[,2],x$Data$Males[,2],x$Data$Unsexed[,2])
           age.tag <- c(rep("Females",length(x$Data$Females[,2])),rep("Males",length(x$Data$Males[,2])),rep("Unsexed",length(x$Data$Unsexed[,2])))
           len     <- c(x$Data$Females[,3],x$Data$Males[,3],x$Data$Unsexed[,3])
           len.tag <- c(rep("Females",length(x$Data$Females[,3])),rep("Males",length(x$Data$Males[,3])),rep("Unsexed",length(x$Data$Unsexed[,3])))
           wgt     <- c(x$Data$Females[,4],x$Data$Males[,4],x$Data$Unsexed[,4])
           wgt.tag <- c(rep("Females",length(x$Data$Females[,4])),rep("Males",length(x$Data$Males[,4])),rep("Unsexed",length(x$Data$Unsexed[,4])))
           if(x$Properties["Start.Month"] == x$Properties["End.Month"])
             {
              par(mfrow=c(2,2),oma=c(2,2,2,1), mar=c(4,4,2,2))
              Ecdf(age,group=age.tag,xlab="",
                   ylab="",col=c("red","blue","black"),label.curves=FALSE)
              mtext(expression("Proportion " <= " age"),2,line=2)
              mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2)
              Ecdf(len,group=len.tag,xlab="",
                   ylab="",col=c("red","blue","black"),label.curves=FALSE)
              mtext(expression("Proportion " <= " length"),2,line=2)
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2)
              plot(x=len[1:length(x$Data$Females[,3])],y=wgt[1:length(x$Data$Females[,4])],
                   xlab="",ylab="",pch=1,
                   xlim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),ylim=c(min(wgt,na.rm=TRUE),max(wgt,na.rm=TRUE)),col="red")
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2)
              mtext(paste("Body Weight (",x$Properties["Units.Weight"],")",sep=""),2,line=2)
              points(x=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                     y=wgt[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                     pch=1,col="blue")
              points(x=len[(length(x$Data$Females[,3])+length(x$Data$Males[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3])+length(x$Data$Unsexed[,3]))],
                     y=wgt[(length(x$Data$Females[,4])+length(x$Data$Males[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4])+length(x$Data$Unsexed[,4]))],
                     pch=1,col="black")
              plot(x=age[1:length(x$Data$Females[,2])],y=len[1:length(x$Data$Females[,3])],
                   xlab="",ylab="",pch=1,
                   xlim=c(min(age,na.rm=TRUE),max(age,na.rm=TRUE)),ylim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),col="red")
              mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2)
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),2,line=2)
              points(x=age[(length(x$Data$Females[,2])+1):(length(x$Data$Females[,2])+length(x$Data$Males[,2]))],
                     y=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                     pch=1,col="blue")
              points(x=age[(length(x$Data$Females[,2])+length(x$Data$Males[,2])+1):(length(x$Data$Females[,2])+length(x$Data$Males[,2])+length(x$Data$Unsexed[,2]))],
                     y=len[(length(x$Data$Females[,3])+length(x$Data$Males[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3])+length(x$Data$Unsexed[,3]))],
                     pch=1,col="black")
              mtext(side=3,text="Unsexed",outer=TRUE,col="black",line=-2)
              mtext(side=3,text="Females",outer=TRUE,col="red",line=-1)
              mtext(side=3,text="Males",outer=TRUE,col="blue",line=0)
             }
           else
             {
              gsi     <- c(100*x$Data$Females[,5]/x$Data$Females[,4],100*x$Data$Males[,5]/x$Data$Males[,4],100*x$Data$Unsexed[,5]/x$Data$Unsexed[,4])
              gsi.tag <- c(rep("Females",length(x$Data$Females[,5])),rep("Males",length(x$Data$Males[,5])),rep("Unsexed",length(x$Data$Unsexed[,5])))
              mon     <- c(as.numeric(format(x$Data$Females$Date, "%m")),as.numeric(format(x$Data$Males$Date, "%m")),as.numeric(format(x$Data$Unsexed$Date, "%m")))
              mon.tag <- c(rep("Females",length(x$Data$Females$Date)),rep("Males",length(x$Data$Males$Date)),rep("Unsexed",length(x$Data$Unsexed$Date)))
              layout(matrix(c(1,2,3,4,5,5),3,2, byrow = TRUE))
              par(oma=c(2,2,2,1), mar=c(4,4,2,2))
              Ecdf(age,group=age.tag,xlab="",
                   ylab="",col=c("red","blue","black"),label.curves=FALSE)
              mtext(expression("Proportion " <= " age"),2,line=2)
              mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2)
              Ecdf(len,group=len.tag,xlab="",
                   ylab="",col=c("red","blue","black"),label.curves=FALSE)
              mtext(expression("Proportion " <= " length"),2,line=2)
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2)
              plot(x=len[1:length(x$Data$Females[,3])],y=wgt[1:length(x$Data$Females[,4])],
                   xlab="",ylab="",pch=1,
                   xlim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),ylim=c(min(wgt,na.rm=TRUE),max(wgt,na.rm=TRUE)),col="red")
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2)
              mtext(paste("Body Weight (",x$Properties["Units.Weight"],")",sep=""),2,line=2)
              points(x=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                     y=wgt[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                     pch=1,col="blue")
              points(x=len[(length(x$Data$Females[,3])+length(x$Data$Males[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3])+length(x$Data$Unsexed[,3]))],
                     y=wgt[(length(x$Data$Females[,4])+length(x$Data$Males[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4])+length(x$Data$Unsexed[,4]))],
                     pch=1,col="black")
              plot(x=age[1:length(x$Data$Females[,2])],y=len[1:length(x$Data$Females[,3])],
                   xlab="",ylab="",pch=1,
                   xlim=c(min(age,na.rm=TRUE),max(age,na.rm=TRUE)),ylim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),col="red")
              mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2)
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),2,line=2)
              points(x=age[(length(x$Data$Females[,2])+1):(length(x$Data$Females[,2])+length(x$Data$Males[,2]))],
                     y=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                     pch=1,col="blue")
              points(x=age[(length(x$Data$Females[,2])+length(x$Data$Males[,2])+1):(length(x$Data$Females[,2])+length(x$Data$Males[,2])+length(x$Data$Unsexed[,2]))],
                     y=len[(length(x$Data$Females[,3])+length(x$Data$Males[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3])+length(x$Data$Unsexed[,3]))],
                     pch=1,col="black")
              plot(1:12,rep(0,12),ylim=c(0,max(gsi,na.rm=TRUE)),ylab="", xlab="",type="n")
              mtext("GSI (%)",2,line=2)
              mtext("Month",1,line=2)
              points(x=mon[1:length(x$Data$Females[,4])],y=gsi[1:length(x$Data$Females[,4])], pch=1, col="red")
              points(x=mon[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                     y=gsi[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                     pch=1, col="blue")
              points(x=mon[(length(x$Data$Females[,4])+length(x$Data$Males[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4])+length(x$Data$Unsexed[,4]))],
                     y=gsi[(length(x$Data$Females[,4])+length(x$Data$Males[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4])+length(x$Data$Unsexed[,4]))],
                     pch=1, col="black")
              mtext(side=3,text="Unsexed",outer=TRUE,col="black",line=-2)
              mtext(side=3,text="Females",outer=TRUE,col="red",line=-1)
              mtext(side=3,text="Males",outer=TRUE,col="blue",line=0)
             }
          }
        else if(!is.na(x$Properties["Units.Weight"]) & !is.na(x$Properties["Units.Gonad"]) & !is.na(x$Properties["Start.Month"]) & !is.na(x$Properties["Units.Liver"]))
          {
           age     <- c(x$Data$Females[,2],x$Data$Males[,2],x$Data$Unsexed[,2])
           age.tag <- c(rep("Females",length(x$Data$Females[,2])),rep("Males",length(x$Data$Males[,2])),rep("Unsexed",length(x$Data$Unsexed[,2])))
           len     <- c(x$Data$Females[,3],x$Data$Males[,3],x$Data$Unsexed[,3])
           len.tag <- c(rep("Females",length(x$Data$Females[,3])),rep("Males",length(x$Data$Males[,3])),rep("Unsexed",length(x$Data$Unsexed[,3])))
           wgt     <- c(x$Data$Females[,4],x$Data$Males[,4],x$Data$Unsexed[,4])
           wgt.tag <- c(rep("Females",length(x$Data$Females[,4])),rep("Males",length(x$Data$Males[,4])),rep("Unsexed",length(x$Data$Unsexed[,4])))
           if(x$Properties["Start.Month"] == x$Properties["End.Month"])
             {
              par(mfrow=c(2,2),oma=c(2,2,2,1), mar=c(4,4,2,2))
              Ecdf(age,group=age.tag,xlab="",
                   ylab="",col=c("red","blue","black"),label.curves=FALSE)
              mtext(expression("Proportion " <= " age"),2,line=2)
              mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2)
              Ecdf(len,group=len.tag,xlab="",
                   ylab="",col=c("red","blue","black"),label.curves=FALSE)
              mtext(expression("Proportion " <= " length"),2,line=2)
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2)
              plot(x=len[1:length(x$Data$Females[,3])],y=wgt[1:length(x$Data$Females[,4])],
                   xlab="",ylab="",pch=1,
                   xlim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),ylim=c(min(wgt,na.rm=TRUE),max(wgt,na.rm=TRUE)),col="red")
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2)
              mtext(paste("Body Weight (",x$Properties["Units.Weight"],")",sep=""),2,line=2)
              points(x=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                     y=wgt[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                     pch=1,col="blue")
              points(x=len[(length(x$Data$Females[,3])+length(x$Data$Males[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3])+length(x$Data$Unsexed[,3]))],
                     y=wgt[(length(x$Data$Females[,4])+length(x$Data$Males[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4])+length(x$Data$Unsexed[,4]))],
                     pch=1,col="black")
              plot(x=age[1:length(x$Data$Females[,2])],y=len[1:length(x$Data$Females[,3])],
                   xlab="",ylab="",pch=1,
                   xlim=c(min(age,na.rm=TRUE),max(age,na.rm=TRUE)),ylim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),col="red")
              mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2)
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),2,line=2)
              points(x=age[(length(x$Data$Females[,2])+1):(length(x$Data$Females[,2])+length(x$Data$Males[,2]))],
                     y=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                     pch=1,col="blue")
              points(x=age[(length(x$Data$Females[,2])+length(x$Data$Males[,2])+1):(length(x$Data$Females[,2])+length(x$Data$Males[,2])+length(x$Data$Unsexed[,2]))],
                     y=len[(length(x$Data$Females[,3])+length(x$Data$Males[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3])+length(x$Data$Unsexed[,3]))],
                     pch=1,col="black")
              mtext(side=3,text="Unsexed",outer=TRUE,col="black",line=-2)
              mtext(side=3,text="Females",outer=TRUE,col="red",line=-1)
              mtext(side=3,text="Males",outer=TRUE,col="blue",line=0)
             }
           else
             {
              gsi     <- c(100*x$Data$Females[,5]/x$Data$Females[,4],100*x$Data$Males[,5]/x$Data$Males[,4],100*x$Data$Unsexed[,5]/x$Data$Unsexed[,4])
              gsi.tag <- c(rep("Females",length(x$Data$Females[,5])),rep("Males",length(x$Data$Males[,5])),rep("Unsexed",length(x$Data$Unsexed[,5])))
              hsi     <- c(100*x$Data$Females[,6]/x$Data$Females[,4],100*x$Data$Males[,6]/x$Data$Males[,4],100*x$Data$Unsexed[,6]/x$Data$Unsexed[,4])
              hsi.tag <- c(rep("Females",length(x$Data$Females[,6])),rep("Males",length(x$Data$Males[,6])),rep("Unsexed",length(x$Data$Unsexed[,6])))
              mon     <- c(as.numeric(format(x$Data$Females$Date, "%m")),as.numeric(format(x$Data$Males$Date, "%m")),as.numeric(format(x$Data$Unsexed$Date, "%m")))
              mon.tag <- c(rep("Females",length(x$Data$Females$Date)),rep("Males",length(x$Data$Males$Date)),rep("Unsexed",length(x$Data$Unsexed$Date)))
              layout(matrix(c(1,2,3,4,5,5,6,6),4,2, byrow = TRUE))
              par(oma=c(2,2,2,1), mar=c(4,4,2,2))
              Ecdf(age,group=age.tag,xlab="",
                   ylab="",col=c("red","blue","black"),label.curves=FALSE)
              mtext(expression("Proportion " <= " age"),2,line=2,cex=0.75)
              mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2,cex=0.75)
              Ecdf(len,group=len.tag,xlab="",
                   ylab="",col=c("red","blue","black"),label.curves=FALSE)
              mtext(expression("Proportion " <= " length"),2,line=2,cex=0.75)
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2,cex=0.75)
              plot(x=len[1:length(x$Data$Females[,3])],y=wgt[1:length(x$Data$Females[,4])],
                   xlab="",ylab="",pch=1,
                   xlim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),ylim=c(min(wgt,na.rm=TRUE),max(wgt,na.rm=TRUE)),col="red")
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),1,line=2,cex=0.75)
              mtext(paste("Body Weight (",x$Properties["Units.Weight"],")",sep=""),2,line=2,cex=0.75)
              points(x=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                     y=wgt[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                     pch=1,col="blue")
              points(x=len[(length(x$Data$Females[,3])+length(x$Data$Males[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3])+length(x$Data$Unsexed[,3]))],
                     y=wgt[(length(x$Data$Females[,4])+length(x$Data$Males[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4])+length(x$Data$Unsexed[,4]))],
                     pch=1,col="black")
              plot(x=age[1:length(x$Data$Females[,2])],y=len[1:length(x$Data$Females[,3])],
                   xlab="",ylab="",pch=1,
                   xlim=c(min(age,na.rm=TRUE),max(age,na.rm=TRUE)),ylim=c(min(len,na.rm=TRUE),max(len,na.rm=TRUE)),col="red")
              mtext(paste("Age (",x$Properties["Units.Age"],")",sep=""),1,line=2,cex=0.75)
              mtext(paste(x$Properties["Length.Type"]," Length (",x$Properties["Units.Length"],")",sep=""),2,line=2,cex=0.75)
              points(x=age[(length(x$Data$Females[,2])+1):(length(x$Data$Females[,2])+length(x$Data$Males[,2]))],
                     y=len[(length(x$Data$Females[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3]))],
                     pch=1,col="blue")
              points(x=age[(length(x$Data$Females[,2])+length(x$Data$Males[,2])+1):(length(x$Data$Females[,2])+length(x$Data$Males[,2])+length(x$Data$Unsexed[,2]))],
                     y=len[(length(x$Data$Females[,3])+length(x$Data$Males[,3])+1):(length(x$Data$Females[,3])+length(x$Data$Males[,3])+length(x$Data$Unsexed[,3]))],
                     pch=1,col="black")
              plot(1:12,rep(0,12),ylim=c(0,max(gsi,na.rm=TRUE)),ylab="", xlab="",type="n")
              mtext("GSI (%)",2,line=2,cex=0.75)
              mtext("Month",1,line=2,cex=0.75)
              points(x=mon[1:length(x$Data$Females[,4])],y=gsi[1:length(x$Data$Females[,4])], pch=19, col="red")
              points(x=mon[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                     y=gsi[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                     pch=1, col="blue")
              points(x=mon[(length(x$Data$Females[,4])+length(x$Data$Males[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4])+length(x$Data$Unsexed[,4]))],
                     y=gsi[(length(x$Data$Females[,4])+length(x$Data$Males[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4])+length(x$Data$Unsexed[,4]))],
                     pch=1, col="black")
              plot(1:12,rep(0,12),ylim=c(0,max(hsi,na.rm=TRUE)),ylab="", xlab="",type="n")
              mtext("HSI (%)",2,line=2,cex=0.75)
              mtext("Month",1,line=2,cex=0.75)
              points(x=mon[1:length(x$Data$Females[,4])],y=hsi[1:length(x$Data$Females[,4])], pch=1, col="red")
              points(x=mon[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                     y=hsi[(length(x$Data$Females[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4]))],
                     pch=1, col="blue")
              points(x=mon[(length(x$Data$Females[,4])+length(x$Data$Males[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4])+length(x$Data$Unsexed[,4]))],
                     y=hsi[(length(x$Data$Females[,4])+length(x$Data$Males[,4])+1):(length(x$Data$Females[,4])+length(x$Data$Males[,4])+length(x$Data$Unsexed[,4]))],
                     pch=1, col="black")
              mtext(side=3,text="Unsexed",outer=TRUE,col="black",line=-2)
              mtext(side=3,text="Females",outer=TRUE,col="red",line=-1)
              mtext(side=3,text="Males",outer=TRUE,col="blue",line=0)
             }
          }

       }
    }
