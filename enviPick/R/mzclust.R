mzclust <-
function(           MSlist,
                    dmzdens=10,
                    ppm=TRUE,
                    drtdens=60,
                    minpeak=4,
                    maxint=1E6,
                    progbar=FALSE,
                    merged=TRUE,
                    from=FALSE,
                    to=FALSE
                     ){

    ############################################################################
    if(!is.logical(ppm)){stop("invalid ppm argument, TRUE or FALSE")}
    if(!is.logical(merged)){stop("invalid merged argument, TRUE or FALSE")}
    if(!is.logical(progbar)){stop("invalid progbar argument, TRUE or FALSE")}
    if(!length(MSlist)==8){stop("This is not an MSlist object")}
    if(minpeak<1){stop("invalid minpeak argument")}
    if(dmzdens<=0 || drtdens<=0){stop("invalid: some argument <= 0")}
    if(!is.loaded("getEIC")){stop(".Call to getEIC failed")}
    if(!MSlist[[1]][[2]]){stop("MSlist invalid. Use mzMLtoMSlist and mzpart to upload and partition .mzML data first.")}
    if(ppm){ppm2=1}else{ppm2=2};
    if(merged){merged2=1}else{merged2=2};
    if(!to){n=length(MSlist[[5]][,1])}else{n=to}
    if(!from){m=1}else{m=from}
    ############################################################################
    # reset if rerun ###########################################################
    MSlist[[6]]<-0;
    MSlist[[7]]<-0;
    MSlist[[8]]<-0;
    MSlist[[4]][[2]][,6]<-rep(0,length(MSlist[[4]][[2]][,4]));
    MSlist[[4]][[2]][,7]<-rep(0,length(MSlist[[4]][[2]][,4]));
    ############################################################################
    # cluster ################################################################
    if(progbar==TRUE){  prog<-winProgressBar("EIC  clustering",min=m,max=n);
                        setWinProgressBar(prog, 0, title = "EIC  clustering", label = NULL);}    
    startat<-c(0);    
    for(k in m:n){ 
      if(progbar==TRUE){setWinProgressBar(prog, k, title = paste("EIC  clustering: ID",k," @",(MSlist[[5]][k,2]-MSlist[[5]][k,1]+1)," measurements",sep=""), label = NULL)}
      if((MSlist[[5]][k,2]-MSlist[[5]][k,1]+1)>1){      
        clusters <- .Call("getEIC",
                          as.numeric(MSlist[[4]][[2]][MSlist[[5]][k,1]:MSlist[[5]][k,2],1]),       # mz
                          as.numeric(MSlist[[4]][[2]][MSlist[[5]][k,1]:MSlist[[5]][k,2],3]),       # RT
                          as.numeric(MSlist[[4]][[2]][MSlist[[5]][k,1]:MSlist[[5]][k,2],2]),       # intens
                          as.integer(order(MSlist[[4]][[2]][MSlist[[5]][k,1]:MSlist[[5]][k,2],2],decreasing=TRUE)),  # intensity order
                          as.integer(order(MSlist[[4]][[2]][MSlist[[5]][k,1]:MSlist[[5]][k,2],3],decreasing=FALSE)),  # RT order                          
                          as.numeric(dmzdens),
                          as.integer(ppm2),
                          as.numeric(drtdens),                                                   
                          as.integer(merged2),
                          PACKAGE="enviPick"
                        )
        clusters[,10]<-clusters[,10]+startat;             
        MSlist[[4]][[2]][(MSlist[[5]][k,1]):(MSlist[[5]][k,2]),6]<-clusters[,10]                       
        MSlist[[4]][[2]][(MSlist[[5]][k,1]):(MSlist[[5]][k,2]),]<-(MSlist[[4]][[2]][(MSlist[[5]][k,1]):(MSlist[[5]][k,2]),][order(clusters[,10],decreasing=FALSE),])
        startat<-c(max(clusters[,10]))
      }else{
        MSlist[[4]][[2]][(MSlist[[5]][k,1]):(MSlist[[5]][k,2]),6]<-startat;   
        startat<-startat+1;   
      }
    }    
    if(progbar==TRUE){close(prog)};
    maxit<-max(MSlist[[4]][[2]][,6]);
    if(maxit>0){
		index <- .Call("indexed",
			as.integer(MSlist[[4]][[2]][,6]),
			as.numeric(MSlist[[4]][[2]][,2]),
			as.integer(minpeak),
			as.numeric(maxint),
			as.integer(maxit),
			PACKAGE="enviPick"
		)
    }
    if(any(index[,2]!=0)){
		index<-index[index[,2]!=0,,drop=FALSE];
		colnames(index)<-c("start_ID","end_ID","number_peaks");
		MSlist[[6]]<-index;
		partID<-.Call("partID",
			as.integer(index),
			as.integer(length(MSlist[[4]][[2]][,6])),
			PACKAGE="enviPick"  
		)
		MSlist[[4]][[2]][,6]<-partID;
    }
    ############################################################################
    MSlist[[1]][[3]]<-TRUE;
    MSlist[[1]][[4]]<-FALSE;
    MSlist[[1]][[5]]<-FALSE;
    MSlist[[3]][[4]]<-length(index[,3]);
    MSlist[[3]][[5]]<-length(MSlist[[4]][[2]][MSlist[[4]][[2]][,6]!=0,6]);
    MSlist[[4]][[2]][,4]<-seq(1,length(MSlist[[4]][[2]][,4]),1);
    ############################################################################
    MSlist[[2]][[2]][14]<-as.character(dmzdens)
    MSlist[[2]][[2]][15]<-as.character(ppm)
    MSlist[[2]][[2]][16]<-as.character(drtdens)
	MSlist[[2]][[2]][17]<-as.character(minpeak)
    MSlist[[2]][[2]][18]<-as.character(maxint)	
    MSlist[[2]][[2]][19]<-as.character(merged)
    MSlist[[2]][[2]][20]<-as.character(from)
    MSlist[[2]][[2]][21]<-as.character(to)	
	############################################################################
    return(MSlist)
}

        