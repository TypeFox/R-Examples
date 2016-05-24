mzpick <-
function(           MSlist,
                    minpeak=4,
                    drtsmall=20,
                    drtfill=10,
                    drttotal=200,
                    recurs=4,
                    weight=2,
                    SB=3,
                    SN=2,
                    minint=1E4,
                    maxint=1E7,
                    ended=2,
                    progbar=FALSE,
                    from=FALSE,
                    to=FALSE
                    ){

    ############################################################################
    # check inputs #############################################################
    if(!length(MSlist)==8){stop("This is not an MSlist object")}
    if(minpeak<=0){stop("minpeak must be >0!")};
    if(drtsmall<=0 || drtsmall<=0){stop("drt must be >0!")};
    if(!is.loaded("pickpeak")){stop(".Call to pickpeak failed!")};
    if(!is.loaded("gapfill")){stop(".Call to gapfill failed!")};
    if(!is.loaded("picklist")){stop(".Call to picklist failed!")};
    if(minint>=maxint){stop("Revise your minint & maxint settings!")}
    if(ended<1){stop("Wrong ended argument!")}
    if(!is.integer(recurs)){recurs<-ceiling(recurs);recurs<-as.integer(recurs)}
    if(!is.integer(ended)){ended<-ceiling(ended);ended<-as.integer(ended)}
    if(!is.integer(minpeak)){minpeak<-ceiling(minpeak);minpeak<-as.integer(minpeak)}
    if(!is.list(MSlist)){stop("Invalid MSlist argument!")}
    if(!MSlist[[1]][[3]]){stop("MSlist invalid. Use mzMLtoMSlist, mzpart and mzclust first.")}
    if(!to){n=length(MSlist[[6]][,1])}else{n=to}
    if(!from){m=1}else{m=from}    
    ############################################################################
    # reset if rerun ###########################################################
    MSlist[[7]]<-0;
    MSlist[[8]]<-0;
    MSlist[[4]][[2]][,7]<-rep(0,length(MSlist[[4]][[2]][,4]));
    ############################################################################
    if(progbar==TRUE){
      prog<-winProgressBar("EIC peak detection",min=m,max=n);
      setWinProgressBar(prog, 0, title = "EIC peak detection", label = NULL);
    }
    startat<-c(0);
    peaknumb<-c(0);
    MSlist[[4]][[2]][,4]<-seq(1,length(MSlist[[4]][[2]][,4]),1); 
    for(i in m:n){    
		if(MSlist[[6]][i,3]>=minpeak){
		  if(progbar==TRUE){setWinProgressBar(prog, i, title = "EIC  peak detection", label = NULL);}
		  # same RT, interpolate ###################################################
		  out1 <- .Call("gapfill",
			as.numeric(MSlist[[4]][[2]][MSlist[[6]][i,1]:MSlist[[6]][i,2],3]),
			as.numeric(MSlist[[4]][[2]][MSlist[[6]][i,1]:MSlist[[6]][i,2],2]),
			as.integer(order(MSlist[[4]][[2]][MSlist[[6]][i,1]:MSlist[[6]][i,2],3],decreasing=FALSE)),
			as.numeric(MSlist[[4]][[2]][MSlist[[6]][i,1]:MSlist[[6]][i,2],1]),
			as.numeric(MSlist[[4]][[2]][MSlist[[6]][i,1]:MSlist[[6]][i,2],4]),
			as.numeric(MSlist[[4]][[1]]),
			as.numeric(drtfill),
			PACKAGE="enviPick"
			)
		  out1<-matrix(out1,ncol=10);
		  colnames(out1)<-c("m/z","intens","RT","index","intens_filt","1pick","pickcrit","baseline","intens_corr","2pick");                       
		  # Filter step ############################################################
		  out1[,5]<-out1[,2];   # yet to be implemented!
		  # 1st peak detection & baseline subtraction & 2nd peak detection ##########  
		  out2 <- .Call("pickpeak",          
			as.numeric(out1),
			as.numeric(drtsmall),
			as.numeric(drttotal),
			as.integer(minpeak),
			as.integer(recurs),
			as.numeric(weight),   # weight
			as.numeric(SB),       # SB 
			as.numeric(SN),       # SN                                                                           
			as.numeric(minint),   # minimum intensity
			as.numeric(maxint),   # maximum intensity threshold
			as.integer(ended),
			as.integer(2),
			PACKAGE="enviPick"
		  )     
		  out2<-matrix(out2,ncol=10);
		  colnames(out2)<-c("m/z","intens","RT","index","intens_filt","1pick","pickcrit","baseline","intens_corr","2pick");
		  # assign final entries ####################################################
		  if(!all(out2[,10]==0)){      
			peaknumb<-peaknumb+max(out2[,10]);
			out2[,10]<-out2[,10]+startat;         
			for(k in 1:length(out2[,10])){
			  if(out2[k,10]!=startat){
				MSlist[[4]][[2]][out2[k,4],7]<-out2[k,10]
			  }
			}        
			startat<-c(max(out2[,10]));
			MSlist[[4]][[2]][MSlist[[6]][i,1]:MSlist[[6]][i,2],]<-
			  MSlist[[4]][[2]][MSlist[[6]][i,1]:MSlist[[6]][i,2],][order(
				MSlist[[4]][[2]][MSlist[[6]][i,1]:MSlist[[6]][i,2],7],decreasing=FALSE),];   
		  }
		}
    }              
    ############################################################################
	maxit<-max(MSlist[[4]][[2]][,7]);
    # generate peak ID table ###################################################
    if(maxit>0){
		if(progbar==TRUE){setWinProgressBar(prog, i, title = "Generate index", label = NULL);}
		index <- .Call("indexed",
			as.integer(MSlist[[4]][[2]][,7]),
			as.numeric(MSlist[[4]][[2]][,2]),
			as.integer(minpeak),
			as.numeric(maxint),
			as.integer(maxit),
			PACKAGE="enviPick"
		)
		if(any(index[,2]!=0)){
			index<-index[index[,2]!=0,,drop=FALSE];
			partID<-.Call("partID",
				as.integer(index),
				as.integer(length(MSlist[[4]][[2]][,7])),
				PACKAGE="enviPick"  
			)
			MSlist[[4]][[2]][,7]<-partID
			colnames(index)<-c("start_ID","end_ID","number_peaks");
			MSlist[[7]]<-index
		}
	}
    ############################################################################
    maxit<-max(MSlist[[4]][[2]][,7]);
	# generate peaklist ########################################################
    if(maxit>0){
 	  peaklist<-matrix(0,ncol=11,nrow=maxit)
      colnames(peaklist)<-c("m/z","var_m/z","max_int","sum_int","RT","minRT","maxRT","part_ID","EIC_ID","peak_ID","Score")
      if(progbar==TRUE){setWinProgressBar(prog, i, title = "Generate peak table", label = NULL);}
      for(i in 1:length(MSlist[[7]][,1])){
        peaklist[i,1]<-mean(MSlist[[4]][[2]][MSlist[[7]][i,1]:MSlist[[7]][i,2],1])
        peaklist[i,2]<-var(MSlist[[4]][[2]][MSlist[[7]][i,1]:MSlist[[7]][i,2],1])
        peaklist[i,3]<-max(MSlist[[4]][[2]][MSlist[[7]][i,1]:MSlist[[7]][i,2],2])
        peaklist[i,4]<-sum(MSlist[[4]][[2]][MSlist[[7]][i,1]:MSlist[[7]][i,2],2])
        peaklist[i,5]<-(MSlist[[4]][[2]][MSlist[[7]][i,1]:MSlist[[7]][i,2],3][MSlist[[4]][[2]][MSlist[[7]][i,1]:MSlist[[7]][i,2],2]==max(MSlist[[4]][[2]][MSlist[[7]][i,1]:MSlist[[7]][i,2],2])[1]])[1]
		peaklist[i,6]<-min(MSlist[[4]][[2]][MSlist[[7]][i,1]:MSlist[[7]][i,2],3])
        peaklist[i,7]<-max(MSlist[[4]][[2]][MSlist[[7]][i,1]:MSlist[[7]][i,2],3])        
        peaklist[i,8:10]<-MSlist[[4]][[2]][MSlist[[7]][i,1],5:7]
        peaklist[i,11]<-0;
      }
      peaklist<-peaklist[order(peaklist[,3],decreasing=TRUE),];
      MSlist[[8]]<-peaklist;
	  MSlist[[3]][[6]]<-length(peaklist[,1])
    }else{
	  MSlist[[3]][[6]]<-"No peaks picked"
	}
    if(progbar==TRUE){close(prog);}
    ############################################################################
    MSlist[[1]][[5]]<-TRUE;
    MSlist[[3]][[7]]<-length(MSlist[[4]][[2]][MSlist[[4]][[2]][,7]!=0,2])  
    ############################################################################
	MSlist[[2]][[2]][22]<-as.character(minpeak)
	MSlist[[2]][[2]][23]<-as.character(drtsmall) 	
	MSlist[[2]][[2]][24]<-as.character(drtfill) 	
	MSlist[[2]][[2]][25]<-as.character(drttotal) 	
	MSlist[[2]][[2]][26]<-as.character(recurs) 	
	MSlist[[2]][[2]][27]<-as.character(weight) 	
	MSlist[[2]][[2]][28]<-as.character(SB) 	
	MSlist[[2]][[2]][29]<-as.character(SN) 	
	MSlist[[2]][[2]][30]<-as.character(minint) 	
	MSlist[[2]][[2]][31]<-as.character(maxint) 	
	MSlist[[2]][[2]][32]<-as.character(ended) 	
	MSlist[[2]][[2]][33]<-as.character(from) 	
	MSlist[[2]][[2]][34]<-as.character(to) 	
    ############################################################################
    return(MSlist);

}



             