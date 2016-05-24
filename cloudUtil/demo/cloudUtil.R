library(cloudUtil)
data(cloudms2)

cloudUtilDemoColorMap<-function(){
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
    return(cm)
}

cloudUtil.ProcessTime<-function(cloudms2){
    EC2_1<-cloudms2[cloudms2$CLOUD=="EC2_1",]
    EC2_2<-cloudms2[cloudms2$CLOUD=="EC2_2",]
    cloudms2.m<-merge(EC2_1, EC2_2, by="id")
    x<-(cloudms2.m$END_PROCESS.x-cloudms2.m$BEGIN_PROCESS.x) / 3600
    y<-(cloudms2.m$END_PROCESS.y-cloudms2.m$BEGIN_PROCESS.y) / 3600

    plot(y~x,
        xlab="EC2_1",
        ylab="EC2_2",
        main="ProcessTime"
    ); 
    abline(xy.lm<-lm(y~x),col="cornflowerblue",lwd=3)
    summary(xy.lm)
}

cloudUtil.CopyInputTime<-function(cloudms2){
    EC2_1<-cloudms2[cloudms2$CLOUD=="EC2_1",]
    EC2_2<-cloudms2[cloudms2$CLOUD=="EC2_2",]
    cloudms2.m<-merge(EC2_1, EC2_2, by="id")
    x<-(cloudms2.m$END_COPYINPUT.x-cloudms2.m$BEGIN_COPYINPUT.x) / 3600
    y<-(cloudms2.m$END_COPYINPUT.y-cloudms2.m$BEGIN_COPYINPUT.y) / 3600

    plot(y~x,
        xlab="EC2_1",
        ylab="EC2_2",
        main="CopyInputTime",
        sub="log scale",
        log='xy'); 


    p.lm<-lm(y~x)
    xxx<-seq(min(x,na.rm=T),max(x, na.rm=T),length=100)
    yyy<-p.lm$coef[2]*xxx+p.lm$coef[1]
    lines(xxx,yyy,
        col="cornflowerblue",
        lwd=3)

    summary(p.lm)
}

cm<-cloudUtilDemoColorMap()

cloudUtil.ProcessTime(cloudms2)
cloudUtil.CopyInputTime(cloudms2)


boxplot((cloudms2$END_PROCESS-cloudms2$BEGIN_PROCESS)/3600~cloudms2$CLOUD, 
    main="process time",
    col=cm,
    ylab="time [hours]")

boxplot(cloudms2$MZXMLFILESIZE*10^-6/(cloudms2$END_COPYINPUT-cloudms2$BEGIN_COPYINPUT)~cloudms2$CLOUD, 
    main="copy input network throughput",
    col=cm,
    ylab="MBytes/s")


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
        plotConcurrent=TRUE)

cloudUtilPlot(begin=cloudms2$BEGIN_PROCESS, 
        end=cloudms2$END_PROCESS, 
        id=cloudms2$id, 
        group=cloudms2$CLOUD, 
        colormap=cm, 
        normalize=TRUE, 
        plotConcurrent=TRUE,
        plotConcurrentMax=TRUE)
