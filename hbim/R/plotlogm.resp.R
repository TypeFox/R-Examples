`plotlogm.resp`<-
function(D,YLAB="Efficacy",YLIM=c(0,1),XLIM=c(-2,2),TITLE=""){
    #par(mfrow=c(1,1))
    plot(XLIM,YLIM,axes=FALSE,xlab="Standardized Geometric Mean Antibody Response",
        ylab=YLAB,type="n",main=TITLE)
    axis(1,at=c(-2,-1,0,1,2),labels=as.character(c(.01,.1,1,10,100)))
    axis(2)
    box()
    if (length(D$col1)==1){
        lines(D$mu,D$out1,lty=1,col=D$col1)
    } else{
    for (j in 1:length(D$col1)){
        lines(D$mu,D$out1[,j],lty=1,col=D$col1[j])
    }
    }
    for (j in 1:length(D$col2)){
        lines(D$mu,D$out2[,j],lty=5,col=D$col2[j])
    }
    for (j in 1:length(D$col3)){
        lines(D$mu,D$out3[,j],lty=3,col=D$col3[j])
    }
}

