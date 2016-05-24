
regime.internal <- function(TS, q=c(0.9, 0.1)) {
    
    doy <- as.factor(TS$hdoy)
    
    hyrstart <- as.numeric(subset(TS, TS$hmonth==1)$month[1])
    
    if (hyrstart != 1) {
        mlabels <- c(month.abb[hyrstart:12], month.abb[1:(hyrstart-1)])
    } else {mlabels <- month.abb}
    
    ### initialize array to be filled
    Qdoy<-array(data=NA, c(max(as.numeric(doy)),6))
    colnames(Qdoy)<- c("MaxQ", "MinQ", "MeanQ", "Q90", "Q10", "Median")
    
    ### calculate stats
    Qdoy[,1]<-tapply(TS$Flow, doy, max, na.rm=TRUE)
    Qdoy[,2]<-tapply(TS$Flow, doy, min, na.rm=TRUE)
    Qdoy[,3]<-tapply(TS$Flow, doy, mean, na.rm=TRUE)
    Qdoy[,4]<-tapply(TS$Flow, doy, stats::quantile, q[1], na.rm=TRUE)
    Qdoy[,5]<-tapply(TS$Flow, doy, stats::quantile, q[2], na.rm=TRUE)
    Qdoy[,6]<-tapply(TS$Flow, doy, stats::median, na.rm=TRUE)
    
    ### set up polygon for inter-quantile shading
    mdoy<-as.numeric(unique(doy))
    xx<-c(1:max(as.numeric(doy)),max(as.numeric(doy)):1)
    yy<-c(Qdoy[,4],Qdoy[max(as.numeric(doy)):1,5])
    
    ### create plot
    graphics::par(mar=c(4,4,0,0))
    yl1=expression(paste("Discharge (m" ^{3}, "/s)"))
    graphics::plot(Qdoy[,1], col="#6BAED6", type="p", pch=19, cex=0.5, xlab="", ylab="",
         xaxt="n")#max
    graphics::title(ylab=yl1, line=2)
    graphics::points(Qdoy[,2], col="#6BAED6", type="p", pch=19, cex=0.5) #min
    graphics::polygon(xx, yy, col="gray", border="#3182BD")
    graphics::points(Qdoy[,3],col="#08519C",type="l",lwd=2) #mean
    graphics::points(Qdoy[,6],col="gray50",type="l",lwd=2)
    mlabels <- c(mlabels[1], mlabels[6], mlabels[12])
    graphics::axis(1, at=c(1, 152, 335), labels=mlabels)
    
}
