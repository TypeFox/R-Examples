rm.sat <-
function(
	peaklist,
	dmz=0.3,
	drt=0.3,
	intrat=0.01,
	spar=0.8,
	corcut=0.8,
	plotit=TRUE
){

    ############################################################################
	if(!is.data.frame(peaklist)){stop("peaklist must be a data.frame")}
	if(length(peaklist[1,])>3){stop("peaklist with > 3 columns not allowed")}
	if(!is.numeric(peaklist[,1]) || !is.numeric(peaklist[,2]) || !is.numeric(peaklist[,3]) ){stop("peaklist columns not numeric")}
    ############################################################################
    int1<-c();
    nget<-c();
    int<-c();
    mass<-c();
    ert<-c();
    getrid=TRUE;
    along<-order(peaklist[,2],decreasing=TRUE);
    getit<-seq(1,length(peaklist[,1]),1)
    done<-rep(TRUE,length(peaklist[,1]));
    for(i in 1:length(along)){
    if(done[along[i]]){
          getit1<-getit[
                        (peaklist[,2]/peaklist[along[i],2])<=intrat   &
                        (abs(peaklist[,3]-peaklist[along[i],3]))<=drt &
                        (abs(peaklist[,1]-peaklist[along[i],1]))<=dmz &
                        done
                       ];
          if(length(getit1)>1){
              if(plotit==TRUE){
                int1<-c(int1,peaklist[along[i],2]);
                nget<-c(nget,length(getit1));
                int<-c(int,peaklist[getit1,2]/peaklist[along[i],2]);
                mass<-c(mass,peaklist[getit1,1]-peaklist[along[i],1]);
                ert<-c(ert,peaklist[getit1,3]-peaklist[along[i],3]);
              }
              dat<-peaklist[getit1,]
              dat1<-peaklist[along[i],]
              dm<-c(dat1[,1]-dat[,1])
              dI<-c(abs(dat[,2]/dat1[,2]))
              dm1<-dm[dm>0]
              dm2<-dm[dm<0]
              dI1<-dI[dm>0]
              dI1<-dI1[order(dm1)]
              dm1<-dm1[order(dm1)]
              dm1<-dm1[dI1<=intrat]
              dI1<-dI1[dI1<=intrat]
              dI2<-dI[dm<0]
              dm2<-abs(dm2)
              dI2<-dI2[order(dm2)]
              dm2<-dm2[order(dm2)]
              dm2<-dm2[dI2<=intrat]
              dI2<-dI2[dI2<=intrat]
              if( length(dm1)>=4 & length(dm2)>=4 ){
              while(
                    abs(dm2[length(dm2)]-dm1[length(dm1)])>(min((dm1[length(dm1)]-dm1[1]),(dm2[length(dm2)]-dm2[1])))
                    &
                    length(dm2)>3
                    &
                    length(dm1)>3
                    ){
                     if(dm2[length(dm2)]>dm1[length(dm1)]){
                         dm2<-dm2[-length(dm2)]
                         dI2<-dI2[-length(dI2)]
                     }else{
                         dm1<-dm1[-length(dm1)]
                         dI1<-dI1[-length(dI1)]
                     }
                  }
              }else{
                if(getrid==TRUE){done[getit1]<-FALSE}
              }
              if( length(dm1)>=4 & length(dm2)>=4 ){
                  s1<-stats::smooth.spline(data.frame(x=c(dm1),y=c(dI1)),spar=spar);
                  s2<-stats::smooth.spline(data.frame(x=c(dm2),y=c(dI2)),spar=spar);
                  i1<-max(dm1[1],dm2[1])
                  i2<-min(dm1[length(dm1)],dm2[length(dm2)])
                  these<-seq(i1,i2,(i2-i1)/20);
                  s1<-predict(s1,these);
                  s2<-predict(s2,these);
                  if(cor(s1$y,s2$y,method="pearson")>=corcut){
                      done[getit1]<-FALSE
                  }
              }else{
                  if(getrid==TRUE){done[getit1]<-FALSE}
              }
          } # if getit1>0
    }
    }
    if(plotit==TRUE & length(mass)>0){
      def.par <- par(no.readonly = TRUE) # save default, for resetting...
      par(mar=c(4,5,1,1));
      plot.new();
      sc<-close.screen();if(sc[1]!=FALSE){for(m in 1:length(sc)){close.screen(sc[m])}};
      split.screen(c(3,1));
        screen(1);
			plot(
				mass,int,
				cex=0.4,pch=19,
				ylab="Intensity ratio to \n main peak",
				xlab="m/z distance from main peak");
			abline(v=0,col="red");abline(h=intrat,col="red");
        screen(2);
			plot(
				ert,int,
				cex=0.4,pch=19,
				ylab="Intensity ratio to \n main peak",
				xlab="RT distance from main peak");
			abline(h=intrat,col="red");
        screen(3);
			plot(
				log10(int1),nget,
				cex=0.5,pch=19,
				xlab="log10( intensity main peak )",
				ylab="Number of peaks \n within tolerances")
        close.screen(all.screens = TRUE);
        par(def.par)#- reset to default
    }
    cat(paste("Fraction of peaklist removed: ",round(length(done[done==FALSE])/length(done),digits=3),"\n",sep=""));
    that<-data.frame(peaklist,done)
    names(that)<-c("mass","intensity","rt","not satellite?")
	############################################################################
    return(that)
}















