plot.mctp <-
function(x,...)
{
nc<-length(x$connames)
text.Ci<-paste(x$input$conf.level*100, "%", "Simultaneous Confidence Intervals")
 Lowerp<-"|"
       plot(x$Analysis$Estimator,1:nc,xlim=c(-1,1), pch=15,axes=FALSE,xlab="",ylab="")
       points(x$Analysis$Lower,1:nc, pch=Lowerp,font=2,cex=2)
              points(x$Analysis$Upper,1:nc, pch=Lowerp,font=2,cex=2)
              abline(v=0, lty=3,lwd=2)
              for (ss in 1:nc){
              polygon(x=c(x$Analysis$Lower[ss],x$Analysis$Upper[ss]),y=c(ss,ss),lwd=2)}
              axis(1, at = seq(-1, 1, 0.1))
              axis(2,at=1:nc,labels=x$connames)
                box()
 title(main=c(text.Ci, paste("Type of Contrast:",x$input$type), paste("Method:", x$AsyMethod)))
}
