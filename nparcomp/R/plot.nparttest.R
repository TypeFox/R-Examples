plot.nparttest <-
function(x,...)
{
nc<-length(x$cmpid)
text.Ci<-paste(x$input$conf.level*100, "%", "Confidence Interval for p")
 Lowerp<-"|"
       plot(x$Analysis$Estimator,1:nc,xlim=c(0,1), pch=15,axes=FALSE,xlab="",ylab="")
       points(x$Analysis$Lower,1:nc, pch=Lowerp,font=2,cex=2)
              points(x$Analysis$Upper,1:nc, pch=Lowerp,font=2,cex=2)
              abline(v=0.5, lty=3,lwd=2)
              for (ss in 1:nc){
              polygon(x=c(x$Analysis$Lower[ss],x$Analysis$Upper[ss]),y=c(ss,ss),lwd=2)}
              axis(1, at = seq(0, 1, 0.1))
              axis(2,at=1:nc,labels=x$cmpid,font=2)
                box()
 title(main=c(text.Ci, paste("Method:", x$AsyMethod) ))
}
