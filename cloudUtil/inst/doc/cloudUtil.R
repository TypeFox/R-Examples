### R code from vignette source 'cloudUtil.Rnw'

###################################################
### code chunk number 1: cloudUtil.Rnw:91-94
###################################################
library(cloudUtil)
data(cloudms2)
cloudms2[sort(sample(nrow(cloudms2),10)),c(1,5,6,15)]


###################################################
### code chunk number 2: cloudUtil.Rnw:107-126
###################################################
hist(cloudms2$END_PREPROCESS - cloudms2$BEGIN_PREPROCESS,100)
## 
boxplot((cloudms2$END_PROCESS-cloudms2$BEGIN_PROCESS)/3600~cloudms2$CLOUD, 
    main="process time",
    ylab="time [hours]")
   
## 
throughput<-cloudms2$MZXMLFILESIZE*10^-6/
(cloudms2$END_COPYINPUT-cloudms2$BEGIN_COPYINPUT)

boxplot(throughput~cloudms2$CLOUD, 
    main="copy input network throughput",
    ylab="MBytes/s")
## 

cloudUtilPlot(begin=cloudms2$BEGIN_PROCESS, 
    end=cloudms2$END_PROCESS, 
    id=cloudms2$id, 
    group=cloudms2$CLOUD)


###################################################
### code chunk number 3: cloudUtil.Rnw:133-166
###################################################
#green
col.amazon<-rgb(0.1,0.8,0.1,alpha=0.2)
col.amazon2<-rgb(0.1,0.8,0.1,alpha=0.2)

#blue
col.fgcz<-rgb(0.1,0.1,0.8,alpha=0.2)
col.fgcz2<-rgb(0.1,0.1,0.5,alpha=0.2)

#red
col.uzh<-rgb(0.8,0.1,0.1,alpha=0.2)
col.uzh2<-rgb(0.5,0.1,0.1,alpha=0.2)

cm<-c(col.amazon, col.amazon2, col.fgcz, col.fgcz2, col.uzh, col.uzh2)

jpeg("cloudms2Fig.jpg", 640, 640)
op<-par(mfrow=c(2,1))
cloudUtilPlot(begin=cloudms2$BEGIN_PROCESS, 
    end=cloudms2$END_PROCESS, 
    id=cloudms2$id, 
    group=cloudms2$CLOUD, 
    colormap=cm, 
    normalize=FALSE, 
    plotConcurrent=TRUE); 

cloudUtilPlot(begin=cloudms2$BEGIN_PROCESS, 
    end=cloudms2$END_PROCESS, 
    id=cloudms2$id, 
    group=cloudms2$CLOUD, 
    colormap=cm, 
    normalize=TRUE, 
    plotConcurrent=TRUE,
    plotConcurrentMax=TRUE)
dev.off()


