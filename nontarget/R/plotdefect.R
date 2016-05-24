plotdefect <-
function(pattern,elements=c("Br")){

    ############################################################################
    # check inputs #############################################################
    # if(length(element)>1){stop("Use only one element or make several plots for several elements.")}
    if(elements[1]!=FALSE){
      for(i in 1:length(elements)){if(any(pattern[[10]]==elements[i])==FALSE){stop("This element was not used for istope pattern screening!")}}
    }
    ############################################################################
    
    ############################################################################
    # calculate mass defects ###################################################
    dmz<-c()
    for(i in 1:length(pattern[[1]][,1])){
          #n<-as.numeric(paste("0.",strsplit(as.character(pattern[[1]][i,1]),".",fixed=TRUE)[[1]][2],sep=""))
          #if(n>0.5){n<-c((1-n)*-1);}
          
          n<-c(pattern[[1]][i,1]-round(pattern[[1]][i,1],digits=0));
          dmz<-c(dmz,n);
    }
    sc<-close.screen();if(sc[1]!=FALSE){for(m in 1:length(sc)){close.screen(sc[m])}};
    plot(pattern[[1]][,1],dmz,pch=19,cex=0.3,xlim=c(100,1000),col="darkgrey",xlab="m/z",ylab="Mass defect");
    ############################################################################
    
    ############################################################################
    # add points to plot #######################################################
    if(elements[1]!=FALSE){
		these<-c()
		for(i in 1:length(elements)){
			these<-c(these,as.character(pattern[[9]][,1][as.character(pattern[[9]][,5])==elements[i]]));
		}
		these<-unique(these)
		for(i in 1:length(pattern[[3]][,1])){
			that<-as.numeric(strsplit(as.character(pattern[[3]][i,2]),",")[[1]]);
			for(j in 1:length(that)){
			  this<-strsplit(as.character(pattern[[1]][that[j],8]),"/")[[1]];
			  if(this[1]!="none"){
				for(n in 1:length(this)){
				  if(any(as.character(this[n])==as.character(these))){
					  points(pattern[[1]][that[j],1],dmz[that[j]],pch=19,col="red",cex=0.5)
				  }
				}
			  }
			}
		}
    }
    ############################################################################

}
