plot.nparttestpaired <-
function(x,...)
{
par(mfrow=c(2,1),oma=c(0,0,0,0))
text.Ci<-paste((x$input$conf.level)*100, "%", "Confidence Interval for p")  
 Lowerp<-"|"
 for (i in 1:2){
       plot(x$Analysis[i,2],1,xlim=c(0,1), pch=15,axes=FALSE,xlab="",ylab="")
       points(x$Analysis[i,1],1, pch=Lowerp,font=2,cex=2)
              points(x$Analysis[i,3],1, pch=Lowerp,font=2,cex=2)
              abline(v=0.5, lty=3,lwd=2)
              polygon(x=c(x$Analysis[i,1],x$Analysis[i,3]),y=c(1,1),lwd=2)
              axis(1, at = seq(0, 1, 0.1))
              axis(2,at=1,labels=x$methodvec[i],font=2)   
              if(i==1){title(main=c(text.Ci, " Method: Brunner-Munzel (BM), Permutation (PERM)" ))}
                box()        
      }   
       
}
