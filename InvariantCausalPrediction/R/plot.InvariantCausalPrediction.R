plot.InvariantCausalPrediction <-
function(x,maxshow=50,col1=NULL,col2=NULL,col3=NULL,mar=c(10,4,3,1),lwd=1,...){
    usedvariables <- x$usedvariables

    if(x$modelReject){
        cat("\n Model has been rejected and there are no valid confidence intervals to plot \n\n")
    }else{
        if( (p <- length(x$colnames))>maxshow){
            cat(paste("\n Number of variables (",p,") exceeds chosen limit (",maxshow,") \n ... only showing first ",maxshow," variables "  ,sep=""))
        }else{
            maxshow <- p
        }
        xvec <- range(unlist(CM <- x$Coeff)) + c(-1,1)*qnorm(1-x$alpha/4)*max(unlist(CS <- x$CoeffVar))
        
        if(is.null(col1)) col1 <- hsv(0.05,0.8,0.8,0.25)
        if(is.null(col2)) col2 <- hsv(0.05,0.8,0.8,0.6)
        if(is.null(col3)) col3 <- hsv(0.62,0.8,0.8,1)
        
        par(mar=mar)
        plot( seq(0.5,maxshow+0.5,length=maxshow), seq(xvec[1],xvec[2],length=maxshow),type="n",xlab="",ylab="",axes=FALSE)
        axis(1, at=1:maxshow,labels= x$colnames[1:maxshow],las=2,cex.axis=if(maxshow<10) 1.5 else 0.8)
        axis(2)
        abline(h=0,lwd=1,col=rgb(0.2,0.2,0.2,0.9))
        noisevec <- runif(max(sapply(CM,length)),-1,1)*0.15
        lwd <- max(0.02,min(1, lwd*sqrt(50/length(noisevec)) ))
        for (pc in 1:maxshow){
            if(pc %in% usedvariables){
                if(length(CM[[pc]])>0){
                    points( (pc + noisevec[1:length(CM[[pc]])]), CM[[pc]], col= col2,pch=20,cex=1)
                    for (ic in 1:length(CM[[pc]])){
                        lines( rep((pc + noisevec[1:length(CM[[pc]])])[ic],2), CM[[pc]][ic] +c(-1,1)*qnorm(1-x$alpha/4)* CS[[pc]][ic], col= col1,lwd=lwd)
                    }
                }
                if( any(x$ConfInt[,pc]==0)) points( pc,0,col= col2,pch=20,cex=1)
                cc <- pc+0.18
                lines( rep(cc,2), ra <- x$ConfInt[,pc],lwd=3,col=col3)
                lines( cc +c(-1,1)*0.1, rep(ra[1],2),lwd=3,col=col3)
                lines( cc +c(-1,1)*0.1, rep(ra[2],2),lwd=3,col=col3)
            }
        }
    }
    cat("\n")

}
