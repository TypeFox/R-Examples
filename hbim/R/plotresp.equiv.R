`plotresp.equiv` <-
function(D,XLIM=c(0,1),YLIM=c(1,100),RLAB="Efficacy of",bounds=XLIM,TITLE=""){
    #par(mfrow=c(1,1))
    plot(XLIM,YLIM,log="y",type="n",xlab=paste(RLAB," antigen 1 alone",sep=""),
       ylab="Equivalent increase in antibody",main=TITLE)

    doline<-function(r1,r2,m1=D$mu,m2=D$mu,LTY=1,COL="black"){
        eab<-equiv.ab(r1,m1,r2,m2,npts=100)
        pick<- eab$equiv.eff2>bounds[1] & eab$equiv.eff2<bounds[2]
        lines(eab$y[pick],10^eab$x[pick],lty=LTY,col=COL)
        #npts<-100
        #e1out<-min(r1)+(max(r1)-min(r1))*(1:npts)/(npts+1)
        #eab<-equiv.increase(m1,r1,m2,r2,e1out)
        #pick<- eab$e2>bounds[1] & eab$e2<bounds[2]
        #lines(eab$e1[pick],eab$equiv.increase[pick],lty=LTY,col=COL)
    }
    if (length(D$col1)==1){
        for (j in 1:length(D$col2)){
            doline(D$out1,D$out2[,j],LTY=5,COL=D$col2[j])
            doline(D$out1,D$out3[,j],LTY=3,COL=D$col3[j])
        } 
    }
    else {
        for (j in 1:length(D$col2)){
            doline(D$out1[,j],D$out2[,j],LTY=5,COL=D$col2[j])
            doline(D$out1[,j],D$out3[,j],LTY=3,COL=D$col3[j])
        } 
    }
}

