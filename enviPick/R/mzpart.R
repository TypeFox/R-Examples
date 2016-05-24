mzpart <-
function(         MSlist,
                  dmzgap=10,       # used for density calculation!
                  drtgap=500,    
                  ppm=TRUE,
                  minpeak=4,
                  peaklimit=2500,
                  cutfrac=0.1,
                  drtsmall=50,        # also used for density calculation!
                  progbar=FALSE,
                  stoppoints=2E5){                  
                  
    ############################################################################
    # check inputs #############################################################
    if(!length(MSlist)==8){stop("This is not an MSlist object")}
    if(minpeak<=1){stop("minpeak must be >1")}
    if(!is.loaded("masspart")){stop(".Call to masspart failed")}
    if(!is.loaded("rtpart")){stop(".Call to rtpart failed")}
    if(!is.loaded("densvec")){stop(".Call to densvec failed")}
    if(cutfrac>=1 || cutfrac<=0){stop("0<cutfrac<<1")}
    if(ppm==TRUE){ppm2=1}else{ppm2=2};
    if(!MSlist[[1]][[1]]){stop("MSlist empty or invalid. Use readMSdata to upload raw .mzML data first.")}
    ############################################################################  
    MSlist[[5]]<-0;
    MSlist[[6]]<-0;
    MSlist[[7]]<-0;
    MSlist[[8]]<-0;
    MSlist[[4]][[2]][,5]<-rep(0,length(MSlist[[4]][[2]][,4]));
    MSlist[[4]][[2]][,6]<-rep(0,length(MSlist[[4]][[2]][,4]));
    MSlist[[4]][[2]][,7]<-rep(0,length(MSlist[[4]][[2]][,4]));
    inlist<-list(0); 
    inlist[[1]]<-c(1,length(MSlist[[4]][[2]][,1]))  
    donemz_in<-TRUE;
    donert_in<-FALSE;
    doneleng_in<-FALSE;
    doneit<-c(1);
    startcount<-inlist[[1]][[2]]; 
    while(any(donemz_in) || any(donert_in) || any(doneleng_in)){
        if(progbar==TRUE){  prog<-winProgressBar(paste("Partition #",doneit),min=0,max=length(inlist));
                            setWinProgressBar(prog, 0, title = paste("Partition #",doneit), label = NULL) }
        donemz_out<-c();
        donert_out<-c();
        doneleng_out<-c();
        outlist<-list(0); 
        counter<-c(1);      
        for(i in 1:length(inlist)){
          if(progbar==TRUE){setWinProgressBar(prog, i, title = paste("Partition #",doneit), label = NULL)}
          if(   donemz_in[i] ||
                donert_in[i] ||
                doneleng_in[i] ){
            if(donemz_in[i]){ # based on dmzgap ###################################               
                MSlist[[4]][[2]][inlist[[i]][1]:inlist[[i]][2],]<-
                  MSlist[[4]][[2]][inlist[[i]][1]:inlist[[i]][2],][order(
                    MSlist[[4]][[2]][inlist[[i]][1]:inlist[[i]][2],1],decreasing=TRUE),]
                part <- .Call("masspart",
                          as.numeric(MSlist[[4]][[2]][inlist[[i]][1]:inlist[[i]][2],1]),
                          as.numeric(dmzgap),
                          as.integer(ppm2),
                          PACKAGE="enviPick"
                        )
                if(max(part)>1){
                  atseq<-c(inlist[[i]][1]:inlist[[i]][2]);
                  from<-1;
                  to<-1;
                  lengit<-abs(inlist[[i]][2]-inlist[[i]][1]+1);
                  while(to<=lengit){
                    while(part[from]==part[to] & to<=lengit){
                      to<-to+1;
                    }
                    if((to-from+1)>=minpeak){             
                      outlist[[counter]]<-c(atseq[from],atseq[(to-1)]); 
                      donemz_out<-c(donemz_out,FALSE);
                      donert_out<-c(donert_out,TRUE);
                      doneleng_out<-c(doneleng_out,((to-from+1)>peaklimit));
                      counter=counter+1;
                    }                    
                    if(doneit==1 & progbar==TRUE){
                      setWinProgressBar(prog, from/lengit, title = "Partition #1", label = NULL);
                    }
                    from<-to;                    
                  }
                }else{
                  if(progbar==TRUE){setWinProgressBar(prog, i, title = paste("Partition #",doneit), label = NULL)}
                  outlist[[counter]]<-inlist[[i]];
                  donemz_out<-c(donemz_out,FALSE);
                  donert_out<-c(donert_out,TRUE);
                  doneleng_out<-c(doneleng_out,(abs(inlist[[i]][2]-inlist[[i]][1]+1)>peaklimit));
                  counter=counter+1;
                }
            }
            if(donert_in[i]){ # based on drtgap ###################################
                MSlist[[4]][[2]][inlist[[i]][1]:inlist[[i]][2],]<-
                  MSlist[[4]][[2]][inlist[[i]][1]:inlist[[i]][2],][order(
                    MSlist[[4]][[2]][inlist[[i]][1]:inlist[[i]][2],3],decreasing=TRUE),]
                part <- .Call("rtpart",
                          as.numeric(MSlist[[4]][[2]][inlist[[i]][1]:inlist[[i]][2],3]),
                          as.numeric(drtgap),
                          PACKAGE="enviPick"
                        )
                if(max(part)>1){
                  atseq<-c(inlist[[i]][1]:inlist[[i]][2]);
                  from<-1;
                  to<-1;
                  lengit<-abs(inlist[[i]][2]-inlist[[i]][1]+1);
                  while(to<=lengit){
                    while(part[from]==part[to] & to<=lengit){
                      to<-to+1;
                    }
                    if((to-from+1)>=minpeak){
                      outlist[[counter]]<-c(atseq[from],atseq[(to-1)]);
                      donemz_out<-c(donemz_out,TRUE);
                      donert_out<-c(donert_out,FALSE);
                      doneleng_out<-c(doneleng_out,((to-from+1)>peaklimit));
                      counter=counter+1;
                    }
                    from<-to;
                  }
                }else{
                  if(progbar==TRUE){setWinProgressBar(prog, i, title = paste("Partition #",doneit), label = NULL)}
                  outlist[[counter]]<-inlist[[i]];
                  donemz_out<-c(donemz_out,FALSE);
                  donert_out<-c(donert_out,FALSE);
                  doneleng_out<-c(doneleng_out,(abs(inlist[[i]][2]-inlist[[i]][1]+1)>peaklimit));
                  counter=counter+1;
                }
            }
            # omit low-density peaks ###########################################
            if( !donert_in[i] &
                !donemz_in[i] &
                doneleng_in[i]){
                 if(progbar==TRUE){setWinProgressBar(prog, i, title = paste("Partition #",doneit,"- omission required for", 
                  length(MSlist[[4]][[2]][inlist[[i]][1]:inlist[[i]][2],3]) , " measurement points"), label = NULL)}
                 if(length(MSlist[[4]][[2]][inlist[[i]][1]:inlist[[i]][2],3])>stoppoints){                  
                    stop("Number measurements remaining in partition > stoppoints; abort. Try using smaller dmzgap or drtgap.");
                 } 
                 MSlist[[4]][[2]][inlist[[i]][1]:inlist[[i]][2],]<-
                  MSlist[[4]][[2]][inlist[[i]][1]:inlist[[i]][2],][order(
                    MSlist[[4]][[2]][inlist[[i]][1]:inlist[[i]][2],1],decreasing=TRUE),]
                 dens <- .Call("densvec",
                          as.numeric(MSlist[[4]][[2]][inlist[[i]][1]:inlist[[i]][2],3]),
                          as.numeric(MSlist[[4]][[2]][inlist[[i]][1]:inlist[[i]][2],1]),
                          as.numeric(dmzgap),
                          as.integer(ppm2),
                          as.numeric(drtsmall),
                          as.numeric(c(1,1)),
                          PACKAGE="enviPick"
                        )
                 MSlist[[4]][[2]][inlist[[i]][1]:inlist[[i]][2],]<-
                 MSlist[[4]][[2]][inlist[[i]][1]:inlist[[i]][2],][order(dens,decreasing=FALSE),];
                 dens<-dens[order(dens,decreasing=FALSE)];
                 dens2<-dens[ceiling(cutfrac*length(dens))];
                 low<-min(seq(inlist[[i]][1],inlist[[i]][2],1)[dens>dens2]);
                 up<-max(seq(inlist[[i]][1],inlist[[i]][2],1)[dens>dens2]);
                 if((up-low+1)>=minpeak){                 
                    inlist[[i]][1]<-low;
                    inlist[[i]][2]<-up; 
                    outlist[[counter]]<-inlist[[i]]                 
                 }
                 donemz_out<-c(donemz_out,TRUE);
                 donert_out<-c(donert_out,FALSE);
                 doneleng_out<-c(doneleng_out,(abs(inlist[[i]][2]-inlist[[i]][1]+1)>peaklimit));
                 counter=counter+1;
            }
          }else{
            if(progbar==TRUE){setWinProgressBar(prog, i, title = paste("Partition #",doneit), label = NULL)}
            outlist[[counter]]<-inlist[[i]] 
            counter=counter+1;
            donemz_out<-c(donemz_out,FALSE);
            donert_out<-c(donert_out,FALSE);
            doneleng_out<-c(doneleng_out,FALSE);
          }
        }
        donemz_in<-donemz_out;
        donert_in<-donert_out;
        doneleng_in<-doneleng_out;
        inlist<-outlist;
        rm(outlist);
        if(progbar==TRUE){close(prog);}
        doneit<-doneit+1;
    }            
    endcount<-c(0);
    count<-c(1);
    for(i in 1:length(inlist)){
      endcount<-c(endcount+(inlist[[i]][2]-inlist[[i]][1]+1));
      MSlist[[4]][[2]][inlist[[i]][1]:inlist[[i]][2],5]<-count;
      count<-count+1;
    }
    # convert inlist to matrix #################################################
    mat<-matrix(nrow=length(inlist),ncol=2,0)
    for(i in 1:length(inlist)){
      mat[i,]<-inlist[[i]]    
    }
    MSlist[[5]]<-mat;
    MSlist[[4]][[2]][,4]<-seq(1,length(MSlist[[4]][[2]][,4]),1);
    ############################################################################
    MSlist[[3]][[1]]<-startcount;
    MSlist[[3]][[2]]<-length(inlist);
    MSlist[[3]][[3]]<-endcount;
    ############################################################################
    MSlist[[2]][[2]][6]<-as.character(dmzgap)
    MSlist[[2]][[2]][7]<-as.character(drtgap)
	MSlist[[2]][[2]][8]<-as.character(ppm)
	MSlist[[2]][[2]][9]<-as.character(minpeak)
	MSlist[[2]][[2]][10]<-as.character(peaklimit)
	MSlist[[2]][[2]][11]<-as.character(cutfrac)
	MSlist[[2]][[2]][12]<-as.character(drtsmall)
	MSlist[[2]][[2]][13]<-as.character(stoppoints)
    ############################################################################	
    MSlist[[1]][[2]]<-TRUE;
    MSlist[[1]][[3]]<-FALSE;
    MSlist[[1]][[4]]<-FALSE;
    MSlist[[1]][[5]]<-FALSE;
    return(MSlist);
}
