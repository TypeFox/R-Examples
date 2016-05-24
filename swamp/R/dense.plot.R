dense.plot <-
function(feature.assoc,lty=1:2,col=1:2,lwd=c(2,2),ylab="",main="",cex.main=1,cex.lab=1,cex.axis=1){
        yylim<-max(c(density(feature.assoc$observed)$y,density(feature.assoc$permuted)$y))
        if(feature.assoc$method=="correlation"|feature.assoc$method=="numeric")     {xmin=-1  
                                                                                      xmax=1}
        if(feature.assoc$method=="t.test"){xmin=min(c(density(feature.assoc$observed)$x,density(feature.assoc$permuted)$x))
                                  xmax=max(c(density(feature.assoc$observed)$x,density(feature.assoc$permuted)$x))}                          
        if(feature.assoc$method=="AUC"|feature.assoc$method=="linear.model")         {xmin=0
                                                                                      xmax=1}
        if(feature.assoc$method=="correlation"|feature.assoc$method=="numeric") {xlab="Correlation"}
        if(feature.assoc$method=="t.test") {xlab="t.test Statistic"}
        if(feature.assoc$method=="AUC") {xlab="Area Under Curve"}
        if(feature.assoc$method=="linear.model") {xlab="R squared"}
        plot(density(feature.assoc$observed),lty=lty[1],col=col[1],lwd=lwd[1],xlim=c(xmin,xmax),ylim=c(0,yylim),xlab=xlab,ylab=ylab,main=main,cex.main=cex.main,cex.lab=cex.lab,cex.axis=cex.axis)
        lines(density(feature.assoc$permuted),lty=lty[2],col=col[2],lwd=lwd[2])
        }

