plotisotopes <-
function(
	input
){
    
	############################################################################
    if(names(input)[1]=="Components"){
      if(length(input[[2]])<2){
        stop("No isotope grouping available!")
      }
    }
    # for input data set ... ###################################################
    if(names(input)[1]=="Patterns"){
      if(input[[2]][8]==TRUE){cat("plot.isotopes for deter=TRUE: unknown isotopes!")}
      ##########################################################################
      def.par <- par(no.readonly = TRUE) # save default, for resetting...
      sc<-close.screen();if(sc[1]!=FALSE){for(m in 1:length(sc)){close.screen(sc[m])}};
      plot.new();
      plot.window(xlim=c(0,(length(input[[9]][,4])+1)),ylim=c(0,max((input[[9]][,3]))));
      axis(2);
      title(ylab="Absolute frequency",xlab="Isotope,z");
      for(i in 1:length(input[[9]][,4])){
          points(i,input[[9]][i,3],type="h",lwd=5);
          points(i,input[[9]][i,4],type="h",lwd=5,col="red");
          mtext(input[[9]][i,2],at=i,side=1,cex=0.5,line=-1);
      };
      rm(i);
      it<-floor(length(input[[11]])/2)
      mtext(input[[9]][1,1],at=1,side=1,cex=0.7,line=1)
      for(i in 2:(length(input[[9]][,4]))){
        if(input[[9]][i,1]!=input[[9]][(i-1),1]){
          mtext(input[[9]][i,1],at=i,side=1,cex=0.7,line=1)
          abline(v=(i-0.5),col="lightgrey")
        }
      }
      rm(i);
      legend((length(input[[9]][,4])+1)/3,max((input[[9]][,3])),
      pch=c(19,19),legend=c("Number of hits per isotope","Number of pattern groups with that hit"),col=c("black","red"),bg="white",adj=0);
      par<-def.par
      return(input[[9]])
      ##########################################################################
    }
    ############################################################################

    ############################################################################    
    # ... or for components ####################################################
    if(names(input)[1]=="Components"){
        get1<-c();
        get2<-c();
        get3<-c();
        for(i in 1:length(input[[1]][,17])){
          if(input[[1]][i,17]!="-"){
            getit<-strsplit(as.character(input[[1]][i,17]),"/")[[1]];
            for(j in 1:length(getit)){
                get1<-c(get1,substr(getit[j],1,(nchar(getit[j])-7)))
                get2<-c(get2,substr(getit[j],(nchar(getit[j])-6),nchar(getit[j])))
                get3<-c(get3,i);
            }  
          }  
        }
        lev1<-as.character(levels(as.factor(get1)));
        lev2<-as.numeric(levels(as.factor(get3)));
        dat<-data.frame(lev1,rep(0,length(lev1)),rep(0,length(lev1)))
        names(dat)<-c("isotope","per component counts","per component counts (only within small m/z toler.)");
        for(i in 1:length(lev2)){
              get7<-as.character(levels(as.factor(get1[get3==lev2[i]])));
              for(j in 1:length(get7)){
               dat[dat[,1]==get7[j],2]<-c(dat[dat[,1]==get7[j],2]+1)
               if(
                  any((get2[get3==lev2[i]][get1[get3==lev2[i]]==get7[j]])=="(small)")
               ){
                  dat[dat[,1]==get7[j],3]<-c(dat[dat[,1]==get7[j],3]+1)
                }
              }
        }
        ########################################################################
        def.par <- par(no.readonly = TRUE) # save default, for resetting...
        sc<-close.screen();if(sc[1]!=FALSE){for(m in 1:length(sc)){close.screen(sc[m])}};
        plot.new();
        plot.window(xlim=c(0,(length(dat[,2])+1)),ylim=c(0,max((dat[,2]))));
        axis(2);
        title(ylab="Absolute frequency",xlab="Isotope,z=all used");
        for(i in 1:length(dat[,2])){
            points(i,dat[i,2],type="h",lwd=5);
            points(i,dat[i,3],type="h",lwd=5,col="red");
        };
        rm(i);
        it<-0;
        mtext(dat[1,1],at=1,side=1,cex=0.7,line=1)
        for(i in 2:(length(dat[,2]))){
            mtext(dat[i,1],at=i,side=1,cex=0.7,line=1)
            abline(v=(i-0.5),col="lightgrey")
        }
        rm(i);
        legend((length(dat[,2])+1)/3,max((dat[,2])),
        pch=c(19,19),legend=c("For all components","Within small m/z tolerance"),col=c("black","red"),bg="white",adj=0);
        par<-def.par
        return(dat)
        ########################################################################
    }
    ############################################################################

}
