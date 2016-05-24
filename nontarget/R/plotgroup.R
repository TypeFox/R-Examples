plotgroup <-function(
	pattern,
	adduct=FALSE,
	groupID,
	massrange=10,
	allmass=TRUE
){

  ##############################################################################
  # check arguments
  if(any(grepl(paste("/",groupID,"/",sep=""),as.character(pattern[[3]][,1]),fixed=TRUE))==FALSE){stop("Invalid group number!\n")}
  if(massrange<=0){stop("set massrange >0!")}
  ##############################################################################

  ##############################################################################
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  #plot.new();
  sc<-close.screen();if(sc[1]!=FALSE){for(m in 1:length(sc)){close.screen(sc[m])}};
  split.screen(c(4,1));
  # retrieve group data set! ###################################################
  dat1<-pattern[[1]][as.numeric(strsplit(as.character(pattern[[3]][grepl(paste("/",groupID,"/",sep=""),as.character(pattern[[3]][,1]),fixed=TRUE),2]),",")[[1]]),]
  dat1<-dat1[order(dat1[,1],decreasing=FALSE),];
  charge<-as.numeric(
	strsplit(
		as.character(pattern[[3]][
			grepl(paste("/",groupID,"/",sep=""),as.character(pattern[[3]][,1]),fixed=TRUE)
		,3])
	,"/")[[1]]
  );                              # charge level
  # find all peaks within parameter ranges + plot:
  screen(1); ###################################################################
  par(mar=c(1,4,2,2));
  dat2<-pattern[[1]][(pattern[[1]][,3]>=as.numeric(dat1[1,3]+pattern[[2]][1]) & pattern[[1]][,3]<=as.numeric(dat1[1,3]+pattern[[2]][2])),]
  dat2<-dat2[dat2[,1]<=(max(dat1[,1])+massrange) & dat2[,1]>=(min(dat1[,1])-massrange),];
  plot(dat2[,1],dat2[,2],ylim=c(0,max(dat2[,2])),type="h",
    main=paste(pattern[[4]][groupID]," / pattern group:",groupID,sep=""),
    xlab="m/z",ylab="Intensity");
  points(dat1[,1],dat1[,2],type="h",col="red",xlab="m/z",ylab="intensity");
  # plot dat 1 only:
  screen(2); ###################################################################
  par(mar=c(2,4,2,2))
  plot(dat1[,1],dat1[,2],ylim=c(0,max(dat1[,2])),type="h",xlab="m/z",ylab="Intensity",col="red");
  # plot relations:
  screen(3); ###################################################################
  par(mar=c(2,4,2,2));
  this1<-c();this2<-c();this3<-c();this4<-c();this5<-c();
  this5b<-(as.numeric(strsplit(as.character(dat1[1,10]),"/")[[1]]))             # charge
  this1b<-(as.numeric(strsplit(as.character(dat1[1,7]),"/")[[1]]));
  this3b<-((strsplit(as.character(dat1[1,8]),"/")[[1]]));
  this4b<-((strsplit(as.character(dat1[1,9]),"/")[[1]]));
  for(j in 1:length(this5b)){
        if(any(charge==this5b[j])){
            this1<-c(this1,this1b[j]);                                          # to
            this2<-c(this2,dat1[1,4]);                                          # from
            this3<-c(this3,this3b[j]);                                          # by what
            this4<-c(this4,this4b[j]);                                          # tolerance
            this5<-c(this5,this5b[j]);
        }
  }
  for(i in 2:length(dat1[,7])){
          this5b<-(as.numeric(strsplit(as.character(dat1[i,10]),"/")[[1]]))
          this1b<-(as.numeric(strsplit(as.character(dat1[i,7]),"/")[[1]]));
          this3b<-((strsplit(as.character(dat1[i,8]),"/")[[1]]));
          this4b<-((strsplit(as.character(dat1[i,9]),"/")[[1]]));
            for(j in 1:length(this5b)){
                if(any(charge==this5b[j])){
                    this1<-c(this1,this1b[j]);
                    this3<-c(this3,this3b[j]);
                    this4<-c(this4,this4b[j]);
                    this2<-c(this2,dat1[i,4]);
                    this5<-c(this5,this5b[j]);
                }
            }
  }
  rm(i);
  plot.window(xlim=c(min(dat1[,1]),max(dat1[,1])),ylim=c(0,length(this1)+1));
  #axis(1)
  count<-c(length(this1));
  for(i in 1:length(this1)){
     if(this4[i]=="large"){
        arrows(dat1[dat1[,4]==this2[i],1],count,dat1[dat1[,4]==this1[i],1],count,col="darkgrey",length=0.05)
      }else{
        arrows(dat1[dat1[,4]==this2[i],1],count,dat1[dat1[,4]==this1[i],1],count,col="darkgreen",length=0.05)
      }
  text(dat1[dat1[,4]==this1[i],1],count,labels=paste(this3[i],"/z=",this5[i],sep=""),pos=4,cex=0.6);
  count<-c(count-1);
  }
  rm(i);
  if(any(this4=="small")){
    mtext("within small mass tolerance",side=1,col="darkgreen",at=min(dat1[,1]),adj=0);
    mtext("within large mass tolerance",side=1,col="darkgrey",at=max(dat1[,1]),adj=1)
  }else{
    mtext("all within large mass tolerance",side=1,col="darkgrey");
  }
  screen(4); ###################################################################
  par(mar=c(2,4,2,2));
  if(length(adduct)>1){
    this<-adduct[[1]][dat1[,4],5]
  }else{
    this<-c("0");
  }
  if(any(this!="0")){
      plot.window(xlim=c(min(dat1[,1]),max(dat1[,1])),ylim=c(0,10));
      leng<-c((max(dat1[,1])-min(dat1[,1]))/(length(dat1[,1]))/6)
      for(i in 1:length(this)){
        this2<-as.numeric(strsplit(as.character(this[i]),"/")[[1]]);
        if(all(this2==0)){
          # y-position:
          if(i>1){if((dat1[i,1]-dat1[(i-1),1])<(leng)){here1<-c(here1-2);}else{here1<-c(10)}}else{here1<-c(10)};
          # x-position
          here2<-dat1[i,1];
          text(here2,here1,paste("(",i,"):"," no adduct",sep=""),col="darkblue",cex=0.7);
        }else{
          # y-position:
          if(i>1){if((dat1[i,1]-dat1[(i-1),1])<(leng)){here1<-c(here1-2);}else{here1<-c(10)}}else{here1<-c(10)};
          # x-position
          here2<-dat1[i,1];
          for(j in 1:length(this2)){
            text(here2,here1,paste("(",i,"):"," adduct group: ",adduct[[3]][this2[j],][,1],sep=""),col="darkblue",cex=0.7);
            this10<-as.numeric(strsplit(as.character(adduct[[3]][this2[j],][,2]),",")[[1]]);
            this11<-strsplit(as.character(adduct[[3]][this2[j],][,3]),"/")[[1]];
            others<-c();for(k in 1:length(this11)){if(this10[k]!=dat1[i,4]){others<-paste(others,this11[k],sep=", ")}};
            others<-substr(others,2,nchar(others));
            text(here2,here1-1.3,paste(this11[this10==dat1[i,4]]," vs. ",others,sep=""),col="darkblue",cex=0.7);
            if(length(this2)>1){segments((here2-leng),here1-2.3,(here2+leng),here1-2.3,col="darkblue")};
            here1<-c(here1-3);
          } # for j
        }   # if/else
      }     # for i
      rm(i);
      mtext("(peak number): adduct group: number  // this adduct vs. further adducts",side=1,col="darkblue");
  }else{
      plot.window(xlim=c(0,1),ylim=c(0,1));
      text(0.5,0.5,"no adducts found",col="darkblue");
  }
  sc<-close.screen();if(sc[1]!=FALSE){for(m in 1:length(sc)){close.screen(sc[m])}};
  par<-def.par;#- reset to default
  if(allmass==TRUE){
    dat2<-dat2[order(dat2[,1]),];return(dat2);
  }else{
    return(dat1);
  };
  
}
