`plotresp.mix` <-
function(D,RLAB="Efficacy of",XYLIM=c(0,1),TITLE=""){
    #par(mfrow=c(1,1))
    plot(XYLIM,XYLIM,type="n",xlab=paste(RLAB," antigen 1 alone",sep=""),
        ylab=paste(RLAB," Mixture",sep=""),main=TITLE)

    doline<-function(r1,r2,m1=D$mu,m2=D$mu,LTY=1,COL="black"){
        eab<-equiv.ab(r1,m1,r2,m2,npts=100)
        lines(eab$equiv.eff1,eab$equiv.eff2,lty=LTY,col=COL)
    }
    lines(XYLIM,XYLIM,lty=1,col="grey")
    if (length(D$col1)==1){
        for (j in 1:length(D$col2)){
            doline(D$out1,D$out2[,j],LTY=5,COL=D$col2[j])
            doline(D$out1,D$out3[,j],LTY=3,COL=D$col3[j])
        }
    }
    else{
        for (j in 1:length(D$col2)){
            doline(D$out1[,j],D$out2[,j],LTY=5,COL=D$col2[j])
            doline(D$out1[,j],D$out3[,j],LTY=3,COL=D$col3[j])
        }
    }
}

