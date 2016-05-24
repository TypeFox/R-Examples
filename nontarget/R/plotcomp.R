plotcomp <-
function(
	comp,
	compoID,
	peakID=FALSE
){

            ####################################################################
            def.par <- par(no.readonly = TRUE) # save default, for resetting...
            sc<-close.screen();if(sc[1]!=FALSE){for(m in 1:length(sc)){close.screen(sc[m])}};
            if(compoID==FALSE & peakID==FALSE){stop("what? compoID? peakID?")}
            if(compoID!=FALSE & peakID!=FALSE){stop("what? compoID or peakID?")}
            if(peakID!=FALSE & compoID==FALSE){
                i<-c(1);
                found<-FALSE;
                while(i < length(comp[[1]][,1])){
                  these1<-as.numeric(strsplit(as.character(comp[[1]][i,3]),",")[[1]]);
                  if(as.character(comp[[1]][i,5])!="-"){
                    these2<-as.numeric(strsplit(as.character(comp[[1]][i,5]),",")[[1]])
                  }else{
                    these2<-0;
                  };                  
                  if(as.character(comp[[1]][i,7])!="-"){
                    these3<-as.numeric(strsplit(as.character(comp[[1]][i,7]),",")[[1]])
                  }else{
                    these3<-0;
                  };
                  if( any(these1==as.numeric(peakID)) ||
                      any(these2==as.numeric(peakID)) || 
                      any(these3==as.numeric(peakID))
                      ){
                    compoID<-as.numeric(comp[[1]][i,1]);
                    found<-TRUE;
                  }
                  i<-c(i+1);
                }
            if(found==FALSE){stop("peakID not part of a component!")};
            if(length(compoID)>1){
                compi<-as.character(compoID[1]);
                for(j in 2:length(compoID)){
                  compi<-paste(compi,",",as.character(compoID[j]),sep="");
                }
                stop(paste("Several component IDs found for given peakID: ",compi,". Select one and use it as compoID instead of peakID! ",sep=""))
            }
            rm(i);
            }
            ####################################################################
            # get: all peaks concerned and subdataset ##########################
            get1a<-as.numeric(strsplit(as.character(comp[[1]][compoID,3]),",")[[1]]);
            if(as.character(comp[[1]][compoID,5])!="-"){
              get1b<-as.numeric(strsplit(as.character(comp[[1]][compoID,5]),",")[[1]]);
            }else{
              get1b<-c()
            }
            if(as.character(comp[[1]][compoID,7])!="-"){
              get1c<-as.numeric(strsplit(as.character(comp[[1]][compoID,7]),",")[[1]]);
            }else{
              get1c<-c()
            }
            get1<-c(get1a,get1b);
            get1<-as.numeric(as.character(levels(as.factor(get1))));
            get3<-c(get1a,get1b,get1c);
            get3<-as.numeric(as.character(levels(as.factor(get3))));
            if(length(get3)<2){
              plot.new()
              text(0.5,0.5,labels=paste("Component ",compoID," only contains one peak, peakID: ", get3," .\n Hence, no plot available!",sep=""),cex=1)
              stop(paste(compoID," only contains one peak, peakID: ", get3," . \nHence, no plot available!",sep=""))
            }
            if(length(comp[[2]])>1){
              dat1<-comp[[2]][get3,];
              ord<-rank(1/dat1[,2]);
              get3<-get3[order(dat1[,1],decreasing=FALSE)];
              ord<-ord[order(dat1[,1],decreasing=FALSE)];
            }else{
              dat1<-FALSE;
            }
            if(length(comp[[3]])>1){
              dat2<-comp[[3]][get3,];
              ord<-rank(1/dat2[,2]);
              get3<-get3[order(dat2[,1],decreasing=FALSE)];
              ord<-ord[order(dat2[,1],decreasing=FALSE)];              
            }else{
              dat2<-FALSE;
            }
            ####################################################################
            # extract isotope pattern relations for all peaks ##################
            if(length(comp[[2]])>1){
                relat1<-matrix(ncol=length(get3),nrow=length(get3),"");
                rownames(relat1)<-get3;colnames(relat1)<-get3;
                for(i in 1:length(get3)){
                     that1<-as.numeric(strsplit(as.character(dat1[dat1[,4]==get3[i],7]),"/")[[1]]);
                     that2<-strsplit(as.character(dat1[dat1[,4]==get3[i],8]),"/")[[1]];
                     that3<-strsplit(as.character(dat1[dat1[,4]==get3[i],10]),"/")[[1]];
                     that4<-strsplit(as.character(dat1[dat1[,4]==get3[i],9]),"/")[[1]];
                     if(that1[1]!=0){ # enter into matrix
                        for(j in 1:length(that1)){
                            relat1[get3==get3[i],get3==that1[j]]<-paste(relat1[get3==get3[i],get3==that1[j]],"/",that2[j],"(",that4[j],")",",z=",that3[j],sep="");
                        };
                     };
                };
                rm(i);
            }
            ####################################################################
            #extract adduct relations for all peaks ############################
            if(length(comp[[3]])>1){
                relat2<-matrix(ncol=length(get3),nrow=length(get3),"");
                rownames(relat2)<-get3;colnames(relat2)<-get3;
                for(i in 1:length(get3)){
                     that1<-as.numeric(strsplit(as.character(dat2[dat2[,4]==get3[i],6]),"/")[[1]]);
                     that2<-strsplit(as.character(dat2[dat2[,4]==get3[i],7]),"//")[[1]];
                     if(that1[1]!=0){ # enter into matrix
                        for(j in 1:length(that1)){
                            relat2[get3==get3[i],get3==that1[j]]<-paste(relat2[get3==get3[i],get3==that1[j]],that2[j],sep="/");
                        };
                     };
                };
                rm(i);
            }
            ####################################################################
            # plot #############################################################
            layout(matrix(ncol=2,nrow=5,c(1,1,1,1,2,1,1,1,1,2)));
            par(mar=c(1,1,1,1));
            plot.new();
            plot.window(xlim=c(-1.1,1.1),ylim=c(-1.1,1.1));
            #title(main=paste("component #",compoID));
            coordx<-rep(0,length(get3));
            coordy<-rep(0,length(get3));
            b=0;
            a=720/(length(get3)+3)/100;
            for(i in 1:length(get3)){
              coordx[i]=sin(b);
              coordy[i]=cos(b);
              b=b+a;
            }
            rm(i);
            text(coordx,coordy,labels=paste(get3,"(",ord,")"),cex=1.3);
            for(i in 1:length(get3)){
              if(any(get3[i]==get1)){
              text(coordx[i],coordy[i],labels=paste(get3[i],"(",ord[i],")"),col="darkgreen",cex=1.3)
            }}
            rm(i);
            coordx<-coordx*0.8;
            coordy<-coordy*0.8;
            if(length(comp[[2]])>1){
              for(i in 1:length(get3)){for(j in 1:length(get3)){
                 if(relat1[i,j]!=""){
                 if(any(get1==get3[i]) & any(get1==get3[j])){
                  arrows(coordx[i],coordy[i],coordx[j],coordy[j],col="blue",length=.1,lwd=2);
                 }else{
                  arrows(coordx[i],coordy[i],coordx[j],coordy[j],col="blue",length=.1,lwd=1);
                 }
              }};rm(i);};
            }
            if(length(comp[[3]])>1){
              for(i in 1:length(get3)){for(j in 1:length(get3)){
                 if(relat2[i,j]!=""){
                 if(any(get1==get3[i]) & any(get1==get3[j])){
                  lines(c(coordx[i],coordx[j]),c(coordy[i],coordy[j]),col="red",lwd=2);
                 }else{
                  lines(c(coordx[i],coordx[j]),c(coordy[i],coordy[j]),col="red",lwd=1);
                 }
              }};rm(i);};
            }
            # point on most intensive peak #####################################
            if(length(comp[[2]])>1){
              dat3<-comp[[2]][get1,];
              that<-dat3[dat3[,2]==max(dat3[,2]),][,4];
              points(coordx[match(that,get3)],coordy[match(that,get3)],pch=21,cex=3);         
			}
            # mark direction ###################################################
            lines(c(-0.25,-0.21),c(1.15,1.0),col="darkgrey");
            arrows(-0.25,1.15,0.1,1.15,length=0.1,col="darkgrey");
            text(0.2,1.15,labels="m/z",col="darkgrey");
            text(-1,1,labels="isotope relations",col="blue");
            text(-1,0.9,labels="adduct relations",col="red");
            # plot peaks #######################################################
            par(mar=c(4,4,1.5,1.5));
            if(length(comp[[2]])>1 & length(comp[[3]])>1){
              mintol<-c(min(dat1[,3],dat2[,3])-comp[[7]][1]);
              maxtol<-c(max(dat1[,3],dat2[,3])+comp[[7]][1]);
              if(comp[[7]][4]==TRUE){
                minmz<-c(min(dat1[,1],dat2[,1])-(comp[[7]][2]*min(dat1[,1],dat2[,1])/1e6));
                maxmz<-c(max(dat1[,1],dat2[,1])+(comp[[7]][2]*max(dat1[,1],dat2[,1])/1e6));
              }else{
                minmz<-c(min(dat1[,1],dat2[,1])-comp[[7]][2]);
                maxmz<-c(max(dat1[,1],dat2[,1])+comp[[7]][2]);              
              }
            }else{
              if(length(comp[[2]])>1){
                  mintol<-c(min(dat1[,3])-comp[[7]][1]);
                  maxtol<-c(max(dat1[,3])+comp[[7]][1]);
                  if(comp[[6]][4]==TRUE){
                    minmz<-c(min(dat1[,1])-(comp[[7]][2]*min(dat1[,1])/1e6));
                    maxmz<-c(max(dat1[,1])+(comp[[7]][2]*max(dat1[,1])/1e6));
                  }else{
                    minmz<-c(min(dat1[,1])-comp[[7]][2]);
                    maxmz<-c(max(dat1[,1])+comp[[7]][2]);              
                  }
              }
              if(length(comp[[3]])>1){
                  mintol<-c(min(dat2[,3])-comp[[7]][1]);
                  maxtol<-c(max(dat2[,3])+comp[[7]][1]);
                  if(comp[[6]][4]==TRUE){
                    minmz<-c(min(dat2[,1])-(comp[[7]][2]*min(dat2[,1])/1e6));
                    maxmz<-c(max(dat2[,1])+(comp[[7]][2]*max(dat2[,1])/1e6));
                  }else{
                    minmz<-c(min(dat2[,1])-comp[[7]][2]);
                    maxmz<-c(max(dat2[,1])+comp[[7]][2]);              
                  }
              }
            }
            if(length(comp[[2]])>1){
                dat4<-comp[[2]][
                             comp[[2]][,3]>=mintol &
                             comp[[2]][,3]<=maxtol &
                             comp[[2]][,1]>=minmz &
                             comp[[2]][,1]<=maxmz
                              ,]
            }else{
                dat4<-comp[[3]][
                             comp[[3]][,3]>=mintol &
                             comp[[3]][,3]<=maxtol &
                             comp[[3]][,1]>=minmz &
                             comp[[3]][,1]<=maxmz
                              ,]
            }
            plot(dat4[,1],dat4[,2],type="h",xlab="m/z",ylab="Intensity",main=paste("component #",compoID),lwd=2,col="lightgrey");
            if(length(comp[[2]])>2){
              points(dat1[,1],dat1[,2],type="h",lwd=1,col="darkgreen");
            }
            if(length(comp[[3]])>2){
              points(dat2[,1],dat2[,2],type="h",lwd=1,col="darkgreen");
            }
            layout(1);
            par<-def.par #- reset to default
            ####################################################################

            ####################################################################
            # generate relational table ########################################
            these1<-c();
            these2<-c();
            these3<-c();
            for(i in 1:length(get3)){for(j in 1:length(get3)){
                if(length(comp[[2]])>1){
                  if(relat1[i,j]!=""){
                    these1<-c(these1,paste(get3[i],"->",get3[j],sep=""));
                    these2<-c(these2,substr(relat1[i,j],2,nchar(relat1[i,j])));
                    these3<-c(these3,(dat1[dat1[,4]==get3[i],2]/dat1[dat1[,4]==get3[j],2]));
                  };
                };
            }};rm(i);
            for(i in 1:length(get3)){for(j in 1:length(get3)){
                if(length(comp[[3]])>1){                
                  if(relat2[i,j]!=""){
                    these1<-c(these1,paste(get3[i],"-",get3[j],sep=""));
                    these2<-c(these2,substr(relat2[i,j],2,nchar(relat2[i,j])));
                    these3<-c(these3,(dat2[dat2[,4]==get3[i],2]/dat2[dat2[,4]==get3[j],2]));
                  };
                };
            }};rm(i);
            relat1<-data.frame(these1,these2,these3);
            names(relat1)<-c("peaks","relation","intensity ratio");
            if(comp[[1]][compoID,6]!="-"){
              if(length(dat1)>1){
                relat<-list(dat1[order(dat1[,1],decreasing=FALSE),c(4,1,2,3,6)],dat4[order(dat4[,1],decreasing=FALSE),
                c(4,1,2,3)],relat1,as.character(comp[[1]][compoID,6]))
              }else{
                relat<-list(dat2[order(dat2[,1],decreasing=FALSE),c(4,1,2,3,6)],dat4[order(dat4[,1],decreasing=FALSE),
                c(4,1,2,3)],relat1,as.character(comp[[1]][compoID,6]))                
              }
            }else{
              if(length(dat1)>1){ 
                relat<-list(dat1[order(dat1[,1],decreasing=FALSE),c(4,1,2,3,6)],dat4[order(dat4[,1],decreasing=FALSE),
                c(4,1,2,3)],relat1,"Not part of a homologue series")
              }else{
                relat<-list(dat2[order(dat2[,1],decreasing=FALSE),c(4,1,2,3,6)],dat4[order(dat4[,1],decreasing=FALSE),
                c(4,1,2,3)],relat1,"Not part of a homologue series")
              }
              
            }
            names(relat)<-c("a","all peaks within range","relations","Part of homologue series:");
            return(relat);
            ####################################################################

}
