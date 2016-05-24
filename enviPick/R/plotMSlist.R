plotMSlist <-
function(MSlist,RTlimit=FALSE,mzlimit=FALSE,shiny=FALSE,ppmbar=8){

    ############################################################################
    if(!length(MSlist)==8){stop("This is not an MSlist object")}
    if(!MSlist[[1]][[1]]){stop("MSlist initialized, but without data")}
    if(!is.loaded("plotit")){stop(".Call to plotit failed!")};
    if(!is.loaded("binRT")){stop(".Call to binRT failed!")};
    if(!is.loaded("binmz")){stop(".Call to binmz failed!")};
    if(mzlimit[1]==FALSE){
      mzlimit<-c(min(MSlist[[4]][[2]][,1]),max(MSlist[[4]][[2]][,1]))  
    }else{if(!length(mzlimit)==2){stop("invalid mzlimit argument")}}
    if(RTlimit[1]==FALSE){
      RTlimit<-c((min(MSlist[[4]][[1]])-(0.1*min(MSlist[[4]][[1]]))),(max(MSlist[[4]][[1]])+(0.1*max(MSlist[[4]][[1]])))) 
    }else{if(!length(RTlimit)==2){stop("invalid RTlimit argument")}}
    ############################################################################    
        
    ############################################################################
    # Define functions #########################################################
    # initialize plotstate list ################################################
    ini_plotstate<-function(MSlist){
      plotstate<-list(0);
      plotstate[[1]]<-TRUE; # drag
      plotstate[[2]]<-FALSE;# zoom+
      plotstate[[3]]<-FALSE;# zoom-
      plotstate[[4]]<-FALSE;# select
      plotstate[[5]]<-FALSE;# partitioned?
      plotstate[[6]]<-FALSE;# clustered?
      plotstate[[7]]<-FALSE;# picked?
      plotstate[[8]]<-TRUE; # raw?
      plotstate[[9]]<-TRUE; # exit?
      plotstate[[10]]<-TRUE; # TIC? FALSE=base
      plotstate[[11]]<-FALSE; # identify?  
      names(plotstate)<-c("drag","zoomout","zoomin","select","part","clustered","picked","raw","exit","TIC","identif")  
      return(plotstate);
    }
    # set plot area ############################################################
    setplot<-function(){
      par(mar=c(0,0,0,0))
      plot.new();
      plot.window(xlim=c(0,10),ylim=c(0,10));
      #rect(0,0,10,10,border="lightgrey");
      rect(2,1,10,4,border="grey");
      rect(2,5,10,7,border="grey");
      rect(2,8,10,10,border="grey");
      text(6,7.3,"RT",cex=0.8,font=2);
      text(1.3,9,"Intensity",cex=0.8,srt=90,font=2);
      text(6,4.3,"m/z",cex=0.8,font=2);
      text(1.3,6,"Intensity",cex=0.8,srt=90,font=2);
      text(6,.3,"m/z",cex=0.8,font=2);
      text(1.3,3,"RT",cex=0.8,srt=90,font=2);
    }
    # set bar area #############################################################
    setbar<-function(plotstate){
      if(!plotstate$drag){rect(0,9.5,1,10,border="darkgrey",lwd=2,col="lightgrey")}else{rect(0,9.5,1,10,border="darkgrey",lwd=2,col="lightgreen")};text(0.5,9.75,"drag");
      if(!plotstate$zoomout){rect(0,9,1,9.5,border="darkgrey",lwd=2,col="lightgrey")}else{rect(0,9,1,9.5,border="darkgrey",lwd=2,col="lightgreen")};text(0.5,9.25,"zoom +");
      if(!plotstate$zoomin){rect(0,8.5,1,9,border="darkgrey",lwd=2,col="lightgrey")}else{rect(0,8.5,1,9,border="darkgrey",lwd=2,col="lightgreen")};text(0.5,8.75,"zoom -");
      if(!plotstate$select){rect(0,8,1,8.5,border="darkgrey",lwd=2,col="lightgrey")}else{rect(0,8,1,8.5,border="darkgrey",lwd=2,col="lightgreen")};text(0.5,8.25,"select");
      rect(0,7.5,1,8,border="darkgrey",lwd=2,col="lightgrey");text(0.5,7.75,"full view");  
      rect(0,7,1,7.5,border="darkgrey",lwd=2,col="grey");text(0.5,7.25,"help"); 
      if(!plotstate$part){rect(0,4.5,1,5,border="darkgrey",lwd=2,col="lightgrey")}else{rect(0,4.5,1,5,border="darkgrey",lwd=2,col="lightgreen")};text(0.5,4.75,"part");
      if(!plotstate$clustered){rect(0,4,1,4.5,border="darkgrey",lwd=2,col="lightgrey")}else{rect(0,4,1,4.5,border="darkgrey",lwd=2,col="lightgreen")};text(0.5,4.25,"EIC");
      if(!plotstate$picked){rect(0,3.5,1,4,border="darkgrey",lwd=2,col="lightgrey")}else{rect(0,3.5,1,4,border="darkgrey",lwd=2,col="lightgreen")};text(0.5,3.75,"peaks");
      if(!plotstate$raw){rect(0,3,1,3.5,border="darkgrey",lwd=2,col="lightgrey")}else{rect(0,3,1,3.5,border="darkgrey",lwd=2,col="lightgreen")};text(0.5,3.25,"none");
      if(!plotstate$identif){rect(0,2.5,1,3,border="darkgrey",lwd=2,col="lightgrey")}else{rect(0,2.5,1,3,border="darkgrey",lwd=2,col="lightgreen")};text(0.5,2.75,"identif");  
      if(shiny==FALSE){rect(0,0,1,0.5,border="darkgrey",lwd=2,col="orange"); text(0.5,0.25,"exit");}
    }
    # task evaluator ###########################################################
    taskeval<-function(plotstate,at1){
      if( shiny==FALSE){if( at1$x>=0 & at1$x<=1 & at1$y>=0 & at1$y<=0.5  ){plotstate$exit<-FALSE;}}  # quit?  
      if( at1$x>=0 & at1$x<=1 & at1$y>9.5 & at1$y<=10 ){plotstate$drag<-TRUE;plotstate$zoomout<-FALSE;plotstate$zoomin<-FALSE;plotstate$select<-FALSE;}
      if( at1$x>=0 & at1$x<=1 & at1$y>9 & at1$y<=9.5  ){plotstate$drag<-FALSE;plotstate$zoomout<-TRUE;plotstate$zoomin<-FALSE;plotstate$select<-FALSE;}
      if( at1$x>=0 & at1$x<=1 & at1$y>8.5 & at1$y<=9  ){plotstate$drag<-FALSE;plotstate$zoomout<-FALSE;plotstate$zoomin<-TRUE;plotstate$select<-FALSE;}
      if( at1$x>=0 & at1$x<=1 & at1$y>8 & at1$y<=8.5  ){plotstate$drag<-FALSE;plotstate$zoomout<-FALSE;plotstate$zoomin<-FALSE;plotstate$select<-TRUE;}
      if( at1$x>=0 & at1$x<=1 & at1$y>4.5 & at1$y<=5  & MSlist[[1]][[2]]){plotstate$part<-TRUE;plotstate$clustered<-FALSE;plotstate$picked<-FALSE;plotstate$raw<-FALSE;}
      if( at1$x>=0 & at1$x<=1 & at1$y>4 & at1$y<=4.5  & MSlist[[1]][[3]]){plotstate$part<-FALSE;plotstate$clustered<-TRUE;plotstate$picked<-FALSE;plotstate$raw<-FALSE;}
      if( at1$x>=0 & at1$x<=1 & at1$y>3.5 & at1$y<=4  & MSlist[[1]][[5]]){plotstate$part<-FALSE;plotstate$clustered<-FALSE;plotstate$picked<-TRUE;plotstate$raw<-FALSE;}
      if( at1$x>=0 & at1$x<=1 & at1$y>3 & at1$y<=3.5  ){plotstate$part<-FALSE;plotstate$clustered<-FALSE;plotstate$picked<-FALSE;plotstate$raw<-TRUE;}
      if( at1$x>=0 & at1$x<=1 & at1$y>2.5 & at1$y<=3  ){if(plotstate$identif){plotstate$identif<-FALSE}else{plotstate$identif<-TRUE}};
      return(plotstate)
    }
    # data plotting ############################################################
    plotdata<-function(MSlist,mzlimit,RTlimit,plotstate,farb,plot1=TRUE,plot2=TRUE,plot3=TRUE){
        ########################################################################
        if(plotstate$part){whatcolor=1}
        if(plotstate$clustered){whatcolor=2}
        if(plotstate$picked){whatcolor=3}
        if(plotstate$raw){whatcolor=4}
        if(!plotstate$select){text(3,8.75,pos=4,labels="Processing, please wait ... ");}
        out <- .Call(  "plotit",
                        as.numeric(RTlimit[1]),
                        as.numeric(RTlimit[2]),
                        as.numeric(mzlimit[1]),
                        as.numeric(mzlimit[2]),
                        as.numeric(MSlist[[4]][[2]][,1]),
                        as.numeric(MSlist[[4]][[2]][,3]),
                        as.numeric(MSlist[[4]][[2]][,2]),
                        as.numeric(MSlist[[4]][[2]][,5]),
                        as.numeric(MSlist[[4]][[2]][,6]),
                        as.numeric(MSlist[[4]][[2]][,7]),
                        as.integer(whatcolor),
                        PACKAGE="enviPick"                                                                                         
        );
        if(length(out)>0 & length(out[,1])<10000){
            maxint<-out[1,5];
            labelRT<-seq(RTlimit[1],RTlimit[2],((RTlimit[2]-RTlimit[1])/3));
            labelRT<-format(labelRT,digits=4)
            labelint<-seq(0,maxint,((maxint)/2))
            labelint<-format(labelint,scientific=TRUE,digits=2)
            labelmz<-seq(mzlimit[1],mzlimit[2],((mzlimit[2]-mzlimit[1])/3));
            labelmz<-format(labelmz,digits=7)
            if(!plotstate$select){rect(3,8.5,9.5,9,col="white",border="white");}
            # upper plot - RT vs. Int ##########################################
            if(plot1){
              plot.window(  xlim=c((RTlimit[1]-(abs(RTlimit[2]-RTlimit[1])*0.25)),RTlimit[2]),
                          ylim=c((0-(maxint*4*1.1)),maxint*1.1));
              if(!plotstate$raw){
                points(out[,3],out[,2],type="h",col=farb[out[,4]]);             
              }else{
                points(out[,3],out[,2],type="h",col="black");
              }
              axis(1,pos=0,at=as.numeric(labelRT),labels=labelRT,col="grey",cex.axis=.7)
              axis(2,pos=RTlimit[1],at=as.numeric(labelint),labels=labelint,col="grey",cex.axis=.7)                             
            }
            # middle plot - m/z vs. Int ########################################
            if(plot2){
              plot.window(  xlim=c((mzlimit[1]-(abs(mzlimit[2]-mzlimit[1])*0.25)),mzlimit[2]),
                          ylim=c((0-(maxint*5/2*1.1)),0+(maxint*1.1*5/2)));
              if(!plotstate$raw){        
                points(out[,1],out[,2],type="h",col=farb[out[,4]]);    
              }else{
                points(out[,1],out[,2],type="h",col="black");            
              }
              axis(1,pos=0,at=as.numeric(labelmz),labels=labelmz,col="grey",cex.axis=.7)                
              axis(2,pos=mzlimit[1],at=as.numeric(labelint),labels=labelint,col="grey",cex.axis=.7)           
              if( ((mzlimit[2]-mzlimit[1])<0.1)){
                lengppm<-(as.numeric(ppmbar)*mean(mzlimit)/1e6)
                if(lengppm<abs(mzlimit[2]-mzlimit[1])){
                  lines(
                    c((mean(mzlimit)-0.5*lengppm),(mean(mzlimit)+0.5*lengppm)),
                    c(0.5*maxint,0.5*maxint),
                    col="blue",lwd=3)
                }
              }
            }
            # lower plot - m/z vs. RT ##########################################
            if(plot3){
              plot.window(  xlim=c((mzlimit[1]-(abs(mzlimit[2]-mzlimit[1])*0.25)),mzlimit[2]),
                            ylim=c((RTlimit[1]-(1/3*(RTlimit[2]-RTlimit[1]))),(RTlimit[2]+(2*(RTlimit[2]-RTlimit[1])))));                  
              if(!plotstate$raw){        
                if(plotstate$identif){
                    text(out[,1],out[,3],labels=as.character(out[,4]-1),cex=0.7,col=farb[out[,4]]);
                }else{
                    points(out[,1],out[,3],pch=19,cex=0.6,col=farb[out[,4]]);
                }     
              }else{
                points(out[,1],out[,3],pch=19,cex=0.6,col="black");          
              }
              if(plotstate$picked & MSlist[[1]][[5]]){
                points(
                  MSlist[[8]][(MSlist[[8]][,1]<=mzlimit[2]&MSlist[[8]][,1]>=mzlimit[1]&MSlist[[8]][,5]<=RTlimit[2]&MSlist[[8]][,5]>=RTlimit[1]),1] ,
                  MSlist[[8]][(MSlist[[8]][,1]<=mzlimit[2]&MSlist[[8]][,1]>=mzlimit[1]&MSlist[[8]][,5]<=RTlimit[2]&MSlist[[8]][,5]>=RTlimit[1]),5],
                  col=farb[MSlist[[8]][(MSlist[[8]][,1]<=mzlimit[2]&MSlist[[8]][,1]>=mzlimit[1]&MSlist[[8]][,5]<=RTlimit[2]&MSlist[[8]][,5]>=RTlimit[1]),8]+1],pch=0,cex=2);
              }                            
              if( ((mzlimit[2]-mzlimit[1])<0.1)){
                lengppm<-(as.numeric(ppmbar)*mean(mzlimit)/1e6)
                if(lengppm<abs(mzlimit[2]-mzlimit[1])){
                  lines(
                    c((mean(mzlimit)-0.5*lengppm),(mean(mzlimit)+0.5*lengppm)),
                    c(mean(RTlimit),mean(RTlimit)),
                    col="blue",lwd=3)
                }
              }
              axis(1,pos=RTlimit[1],at=as.numeric(labelmz),labels=labelmz,col="grey",cex.axis=.7) 
              axis(2,pos=mzlimit[1],at=as.numeric(labelRT),labels=labelRT,col="grey",cex.axis=.7) 
            }
            ####################################################################
            plot.window(xlim=c(0,10),ylim=c(0,10)); 
        }else{
          if(length(out[,1])>=10000){ # bin data on RT and mz ##################
            if(plotstate$raw){out[,4]<-c(out[,4]+1)} # reset!
            if(plotstate$TIC){what=1}else{what=2}
            RTbin <- .Call (  "binRT",  
                              as.numeric(out[order(out[,3],decreasing=FALSE),3]),
                              as.numeric(out[order(out[,3],decreasing=FALSE),2]),
                              as.numeric(MSlist[[4]][[1]]),
                              as.numeric(out[order(out[,3],decreasing=FALSE),4]), # to select on color
                              as.integer(what),
                              PACKAGE="enviPick"
                            )
            binmzs<-seq(min(out[,1]),max(out[,1]),(max(out[,1])-min(out[,1]))/1000)
            mzbin <- .Call (  "binmz",                
                              as.numeric(out[order(out[,1],decreasing=FALSE),1]),
                              as.numeric(out[order(out[,1],decreasing=FALSE),2]),
                              as.numeric(binmzs),
                              as.numeric(out[order(out[,1],decreasing=FALSE),4]), # to select on color  
                              PACKAGE="enviPick"                         
                           )  
            maxint1<-max(RTbin[,1]);
            maxint2<-max(mzbin[,1]);
            labelRT<-seq(RTlimit[1],RTlimit[2],((RTlimit[2]-RTlimit[1])/3));
            labelRT<-format(labelRT,digits=4)
            labelint1<-seq(0,maxint1,((maxint1)/2))
            labelint1<-format(labelint1,scientific=TRUE,digits=2)
            labelint2<-seq(0,maxint2,((maxint2)/2))
            labelint2<-format(labelint2,scientific=TRUE,digits=2)        
            labelmz<-seq(mzlimit[1],mzlimit[2],((mzlimit[2]-mzlimit[1])/3));
            labelmz<-format(labelmz,digits=7)
            if(!plotstate$select){rect(3,8.5,9.5,9,col="white",border="white");}
            # upper plot - RT vs. Int ##########################################                       
            if(plot1){
              plot.window(  xlim=c((RTlimit[1]-(abs(RTlimit[2]-RTlimit[1])*0.25)),RTlimit[2]),
                            ylim=c((0-(maxint1*4*1.1)),maxint1*1.1));
              points(RTbin[,2],RTbin[,1],type="h",col="black");             
              axis(1,pos=0,at=as.numeric(labelRT),labels=labelRT,col="grey",cex.axis=.7)
              axis(2,pos=RTlimit[1],at=as.numeric(labelint1),labels=labelint1,col="grey",cex.axis=.7)
              plot.window(xlim=c(0,10),ylim=c(0,10));
              text(9.7,9.7,labels="TIC",cex=0.9,col="darkgreen")                             
              }
            # middle plot - m/z vs. Int ########################################
            if(plot2){        
              plot.window(  xlim=c((mzlimit[1]-(abs(mzlimit[2]-mzlimit[1])*0.25)),mzlimit[2]),
                            ylim=c((0-(maxint2*5/2*1.1)),0+(maxint2*1.1*5/2)));
              points(mzbin[,2],mzbin[,1],type="h",col="black");    
              axis(1,pos=0,at=as.numeric(labelmz),labels=labelmz,col="grey",cex.axis=.7)                
              axis(2,pos=mzlimit[1],at=as.numeric(labelint2),labels=labelint2,col="grey",cex.axis=.7)           
            }
            # lower plot - m/z vs. RT ##########################################
            if(plot3){        
              plot.window(xlim=c(0,10),ylim=c(0,10)); 
              text(6,3.3,labels="Too much data points for this panel. Use zoom+ to resolve more details.",cex=0.9,col="darkgreen")
              plot.window(  xlim=c((mzlimit[1]-(abs(mzlimit[2]-mzlimit[1])*0.25)),mzlimit[2]),
                            ylim=c((RTlimit[1]-(1/3*(RTlimit[2]-RTlimit[1]))),(RTlimit[2]+(2*(RTlimit[2]-RTlimit[1])))));                  
              axis(1,pos=RTlimit[1],at=as.numeric(labelmz),labels=labelmz,col="grey",cex.axis=.7) 
              axis(2,pos=mzlimit[1],at=as.numeric(labelRT),labels=labelRT,col="grey",cex.axis=.7) 
            }
          }
        }
    }
    # reset plot limits - within plot operations ###############################
    setlims<-function(mzlimit,RTlimit,at1){
          # lower plot #########################################################
        if(at1$x>=2 & at1$x<=10 & at1$y>=1 & at1$y<=10){
          points(at1$x,at1$y,col="red",pch=19,cex=1.5);
          at2<-locator(n=1);
          if((at2$x>=2 & at2$x<=10 & at2$y>=1 & at2$y<=4)&(at1$x>=2 & at1$x<=10 & at1$y>=1 & at1$y<=4)){
            if(plotstate$drag){ # drag?
              lines(c(at1$x,at2$x),c(at1$y,at2$y),col="red",lwd=2);
              mzlimit<-c(mzlimit-((at2$x-at1$x)/8)*(mzlimit[2]-mzlimit[1]));
              RTlimit<-c(RTlimit-((at2$y-at1$y)/3)*(RTlimit[2]-RTlimit[1]));
              setplot();
              setbar(plotstate);
              plotdata(MSlist,mzlimit,RTlimit,plotstate,farb,plot1=TRUE,plot2=TRUE,plot3=TRUE);              
            }
            if(plotstate$zoomout){ # zoom out?
              rect(at1$x,at1$y,at2$x,at2$y,border="red",lwd=2);          
              delmz<-c(mzlimit[2]-mzlimit[1]);
              mzlimit[2]<-(mzlimit[1]+(((at2$x-2)/8)*delmz));
              mzlimit[1]<-(mzlimit[1]+(((at1$x-2)/8)*delmz));
              mzlimit<-c(min(mzlimit),max(mzlimit));
              delRT<-c(RTlimit[2]-RTlimit[1]);
              RTlimit[2]<-(RTlimit[1]+(((at2$y-1)/3)*delRT));
              RTlimit[1]<-(RTlimit[1]+(((at1$y-1)/3)*delRT));    
              RTlimit<-c(min(RTlimit),max(RTlimit));                      
              setplot();
              setbar(plotstate);
              plotdata(MSlist,mzlimit,RTlimit,plotstate,farb,plot1=TRUE,plot2=TRUE,plot3=TRUE);    
            }
            if(plotstate$zoomin){ # zoom in?      
              rect(at1$x,at1$y,at2$x,at2$y,border="blue",lwd=2);          
              delmz<-c(mzlimit[2]-mzlimit[1]);
              mzlimit[2]<-(mzlimit[2]+(((at2$x-2)/8)*delmz));
              mzlimit[1]<-(mzlimit[1]-(((at1$x-2)/8)*delmz));
              mzlimit<-c(min(mzlimit),max(mzlimit));
              delRT<-c(RTlimit[2]-RTlimit[1]);
              RTlimit[2]<-(RTlimit[2]+(((at2$y-1)/3)*delRT));
              RTlimit[1]<-(RTlimit[1]-(((at1$y-1)/3)*delRT));    
              RTlimit<-c(min(RTlimit),max(RTlimit));                      
              setplot();
              setbar(plotstate);
              plotdata(MSlist,mzlimit,RTlimit,plotstate,farb,plot1=TRUE,plot2=TRUE,plot3=TRUE);    
           }
           if(plotstate$select){ # select
              rect(at1$x,at1$y,at2$x,at2$y,border="black",lwd=1);
              mzlimit2<-mzlimit;
              delmz2<-c(mzlimit2[2]-mzlimit2[1]);
              mzlimit2[2]<-(mzlimit2[1]+(((at2$x-2)/8)*delmz2));
              mzlimit2[1]<-(mzlimit2[1]+(((at1$x-2)/8)*delmz2));
              mzlimit2<-c(min(mzlimit2),max(mzlimit2));          
              RTlimit2<-RTlimit;
              delRT2<-c(RTlimit2[2]-RTlimit2[1]);
              RTlimit2[2]<-(RTlimit2[1]+(((at2$y-1)/3)*delRT2));
              RTlimit2[1]<-(RTlimit2[1]+(((at1$y-1)/3)*delRT2));    
              RTlimit2<-c(min(RTlimit2),max(RTlimit2));                      
              setplot();
              setbar(plotstate);
              plotdata(MSlist,mzlimit2,RTlimit2,plotstate,farb,plot1=TRUE,plot2=TRUE,plot3=FALSE);
              plotdata(MSlist,mzlimit,RTlimit,plotstate,farb,plot1=FALSE,plot2=FALSE,plot3=TRUE);
              rect(at1$x,at1$y,at2$x,at2$y,border="black",lwd=1);
           }
          }
          # middle plot ########################################################
          if((at2$x>=2 & at2$x<=10 & at2$y>=5 & at2$y<=7)&(at1$x>=2 & at1$x<=10 & at1$y>=5 & at1$y<=7)){
            if(plotstate$drag){ # drag?
              lines(c(at1$x,at2$x),c(at1$y,at1$y),col="red",lwd=2); 
              mzlimit<-c(mzlimit-((at2$x-at1$x)/8)*(mzlimit[2]-mzlimit[1]));
            }
            if(plotstate$zoomout){ # zoom out?
              lines(c(at1$x,at2$x),c(at1$y,at1$y),col="red",lwd=2);        
              delmz<-c(mzlimit[2]-mzlimit[1]);
              mzlimit[2]<-(mzlimit[1]+(((at2$x-2)/8)*delmz));
              mzlimit[1]<-(mzlimit[1]+(((at1$x-2)/8)*delmz));
              mzlimit<-c(min(mzlimit),max(mzlimit));                    
            }
            if(plotstate$zoomin){ # zoom in?      
              lines(c(at1$x,at2$x),c(at1$y,at1$y),col="blue",lwd=2);         
              delmz<-c(mzlimit[2]-mzlimit[1]);
              mzlimit[2]<-(mzlimit[2]+(((at2$x-2)/8)*delmz));
              mzlimit[1]<-(mzlimit[1]-(((at1$x-2)/8)*delmz));
              mzlimit<-c(min(mzlimit),max(mzlimit));                     
           }
           setplot();
           setbar(plotstate);
           plotdata(MSlist,mzlimit,RTlimit,plotstate,farb,plot1=TRUE,plot2=TRUE,plot3=TRUE);  
          }
          # upper plot ######################################################### 
          if((at2$x>=2 & at2$x<=10 & at2$y>=8 & at2$y<=10)&(at1$x>=2 & at1$x<=10 & at1$y>=8 & at1$y<=10)){
            if(plotstate$drag){ # drag?
              lines(c(at1$x,at2$x),c(at1$y,at1$y),col="red",lwd=2);
              RTlimit<-c(RTlimit-((at2$x-at1$x)/8)*(RTlimit[2]-RTlimit[1]));
            }
            if(plotstate$zoomout){ # zoom out?
              lines(c(at1$x,at2$x),c(at1$y,at1$y),col="red",lwd=2);         
              delRT<-c(RTlimit[2]-RTlimit[1]);
              RTlimit[2]<-(RTlimit[1]+(((at2$x-2)/8)*delRT));
              RTlimit[1]<-(RTlimit[1]+(((at1$x-2)/8)*delRT));    
              RTlimit<-c(min(RTlimit),max(RTlimit));                      
            }
            if(plotstate$zoomin){ # zoom in?      
              lines(c(at1$x,at2$x),c(at1$y,at1$y),col="blue",lwd=2);          
              delRT<-c(RTlimit[2]-RTlimit[1]);
              RTlimit[2]<-(RTlimit[2]+(((at2$x-2)/8)*delRT));
              RTlimit[1]<-(RTlimit[1]-(((at1$x-2)/8)*delRT));    
              RTlimit<-c(min(RTlimit),max(RTlimit));                      
           }
           setplot();
           setbar(plotstate);
           plotdata(MSlist,mzlimit,RTlimit,plotstate,farb,plot1=TRUE,plot2=TRUE,plot3=TRUE);
         }
      }
      listit<-list(0)
      listit[[1]]<-RTlimit;
      listit[[2]]<-mzlimit;
      return(listit);
    }
    runhelp<-function(){    
      par(mar=c(0,0,0,0))
      plot.new();
      plot.window(xlim=c(0,10),ylim=c(0,10));
      setbar(plotstate)
      rect(0,7,1,7.5,border="red",lwd=2,col="grey");text(0.5,7.25,"help",col="red");
      rect(0,8.5,10,10,border="red",lwd=2);
      text(1,9.25,pos=4,labels=" Zoom and drag within one panel by two-point mouse selections. \n This will update the other panels; panels are always simultaneously active. \n Just like zoom in, extent of zoom out depends on selection scale.",cex=0.8,col="red");     
      rect(0,8,10,8.5,border="red",lwd=2);
      text(1,8.25,pos=4,labels="Selection in lower plot panel (then static) will be shown in the top and middle panels",cex=0.8,col="red");  
      rect(0,7.5,10,8,border="red",lwd=2);
      text(1,7.75,pos=4,labels="Shows all data: TIC (top panel) and highest intensity within m/z bins (middle panel).",cex=0.8,col="red");        
      rect(0,3,10,5,border="red",lwd=2);
      text(1,4,pos=4,labels="All panels: adjust plot color to distinguish either individual partitions, cluster (EICs) or peak \nmeasurements by selecting [part], [EIC] or [peaks], respectively. \nNote that peaks are nested in cluster are nested in partitions. \nThe picked [peaks] are shown as squares.\n[none] depicts the raw data." ,cex=0.8,col="red");
      rect(0,2.5,10,3,border="red",lwd=2);
      text(1,2.75,pos=4,labels="For lower panel: change between scatterplot or IDs of what is set: [part], [EIC] or [peaks].",cex=0.8,col="red");     
      rect(0,0,10,0.5,border="red",lwd=2);
      text(1,0.25,pos=4,labels="Close plot tool and return to R console",cex=0.8,col="red");         
    }
    ############################################################################

    ############################################################################
    # run plotting loop ########################################################
    farb<-colors()[c(157,25,33,53,97,89,81,83,84,69,73,70,71,142,139,556,100,101,102,555,107,417,401,469,490,554,655,657,553)]
    #farb<-c("red","green","blue","orange")
    farb<-sample(farb,max(MSlist[[4]][[2]][,c(5,6,7)]),replace=TRUE)
    farb[1]<-"lightgrey"
    plotstate<-ini_plotstate(MSlist)
    setplot()
    setbar(plotstate)
    plotdata(MSlist,mzlimit,RTlimit,plotstate,farb,plot1=TRUE,plot2=TRUE,plot3=TRUE)                               
    options(locatorBell=FALSE);
    while(plotstate$exit){
      plot.window(xlim=c(0,10),ylim=c(0,10));
      at1<-locator(n=1);
      # (1) refresh settings, and do sth ... 
      if(at1$x<=1){
        plotstate<-taskeval(plotstate,at1)
        setbar(plotstate)
        if(at1$x>=0 & at1$y>2.5 & at1$y<=5){ # recolor data
            setplot()
            setbar(plotstate)
            plotdata(MSlist,mzlimit,RTlimit,plotstate,farb,plot1=TRUE,plot2=TRUE,plot3=TRUE);
            next;
        }
        if(at1$x>=0 & at1$y>7.5 & at1$y<=8){ # full zoom-out
            RTlimit<-c(min(MSlist[[4]][[1]]),max(MSlist[[4]][[1]]))
            mzlimit<-c(min(MSlist[[4]][[2]][,1]),max(MSlist[[4]][[2]][,1]))
            setplot()
            setbar(plotstate)
            plotdata(MSlist,mzlimit,RTlimit,plotstate,farb,plot1=TRUE,plot2=TRUE,plot3=TRUE);
            next;    
        }
        if(at1$x>=0 & at1$y>7 & at1$y<=7.5){
          runhelp();
          at1<-locator(n=1);
          setplot()
          setbar(plotstate)
          plotdata(MSlist,mzlimit,RTlimit,plotstate,farb,plot1=TRUE,plot2=TRUE,plot3=TRUE);          
        }
      # (2) ... or evaluate plot selections for replotting
      }else{
        listit<-setlims(mzlimit,RTlimit,at1)
        RTlimit<-listit[[1]];
        mzlimit<-listit[[2]];
      }
    };
    graphics.off();
    options(locatorBell=TRUE);
    return("done")
    ############################################################################

}
