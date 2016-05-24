readMSdata <-
function( 
			filepath.mzXML,
			MSlevel=c(1),
			progbar=FALSE,
			minRT=FALSE,
			maxRT=FALSE,
			minmz=FALSE,
			maxmz=FALSE
		){
	
		##########################################################################
		# checks #################################################################
		if(!is.character(file.path(filepath.mzXML))){stop("filepath.mzML not a character string")}
		if(!file.exists(file.path(filepath.mzXML))){stop("invalid filepath.mzML")}
		if(minmz==FALSE){
			min_mz<-0
		}else{
			min_mz<-minmz
		}
		if(maxmz==FALSE){
			max_mz<-Inf
		}else{
			max_mz<-maxmz
		}
		##########################################################################
	  
		##########################################################################
		# local function definition ##############################################
		genMSlist<-function(){
			# define MSlist list #################################################
			MSlist<-list(0);
			MSlist[[2]]<-0;MSlist[[3]]<-0;MSlist[[4]]<-0;MSlist[[5]]<-0;MSlist[[6]]<-0;MSlist[[7]]<-0;MSlist[[8]]<-0;
			names(MSlist)<-c("State","Parameters","Results","Scans","Partition_index","EIC_index","Peak_index","Peaklist")
			stage<-data.frame(FALSE,FALSE,FALSE,FALSE,FALSE);
			names(stage)<-c("Raw?","Partitioned?","Clustered?","Filtered?","Picked?");
			MSlist[[1]]<-stage;
			p<-c(   
			"agglom_dmzgap",#1
			"agglom_ppm",
			"agglom_drtgap",
			"agglom_minpeak",
			"agglom_maxint",#5
			"part_dmzgap",#6
			"part_drtgap",
			"part_ppm",
			"part_minpeak",
			"part_peaklimit",
			"part_cutfrac",
			"part_drtsmall",
			"part_stoppoints",#13
			"clust_dmzdens",#14
			"clust_ppm",		  
			"clust_drtdens",
			"clust_minpeak",
			"clust_maxint",			
			"clust_merged",
			"clust_from",
			"clust_to",#21
			"pick_minpeak",#22
			"pick_drtsmall",
			"pick_drtfill",
			"pick_drttotal",
			"pick_recurs",
			"pick_weight",
			"pick_SB",
			"pick_SN",
			"pick_minint",
			"pick_maxint",
			"pick_ended",			
			"pick_from",
			"pick_to"#34	  
			);   
			parameters<-data.frame(p,rep("0",length(p)),stringsAsFactors=FALSE);
			names(parameters)<-c("parameters","value")
			MSlist[[2]]<-parameters;
			######################################################################
			return(MSlist)
		}
		##########################################################################
      
		##########################################################################
		mz1<-readMzXmlFile(file.path(filepath.mzXML),removeMetaData = FALSE,verbose = FALSE)
		peaknumb<-0
		RT<-c()
		done_cen<-FALSE
		for(i in 1:length(mz1)){
			if((any(mz1[[i]]$metaData$msLevel==MSlevel)) & (mz1[[i]]$metaData$peaksCount>0)){
				if((minRT!=FALSE)&(mz1[[i]]$metaData$retentionTime<minRT)) next
				if((maxRT!=FALSE)&(mz1[[i]]$metaData$retentionTime>maxRT)) next
				if(any(names(mz1[[i]]$metaData)=="centroided")){
					if(mz1[[i]]$metaData$centroided!=1){
						stop("\nYour .mzXML-file has not been centroided.\n")
					}
				}else{
					if(!done_cen){
						cat("\nYou have ensured your data is centroided ...\n");
						done_cen<-TRUE
					}
				}
				RT<-c(RT,mz1[[i]]$metaData$retentionTime)	
				peaknumb<-c(peaknumb + sum(mz1[[i]][[1]]$mass>=min_mz & mz1[[i]][[1]]$mass<=max_mz,na.rm=TRUE))
			}
		}
		if(peaknumb==0){
			stop("\nWith this file & settings: no peaks available.\n")
		}
		scans<-list(0)
		scans[[1]]<-RT      
        getpeaks<-matrix(nrow=peaknumb,ncol=7,0)
		from<-1;
		if(progbar==TRUE){  
			prog<-winProgressBar("Extract scans",min=0,max=length(mz1));
            setWinProgressBar(prog, 0, title = "Extract scans", label = NULL);
		}
		for(i in 1:length(mz1)){
			if(progbar==TRUE){setWinProgressBar(prog, i, title = "Extract scans", label = NULL);}
			if(any(mz1[[i]]$metaData$msLevel==MSlevel) & (mz1[[i]]$metaData$peaksCount>0)){
				if((minRT!=FALSE)&(mz1[[i]]$metaData$retentionTime<minRT)) next
				if((maxRT!=FALSE)&(mz1[[i]]$metaData$retentionTime>maxRT)) next
				to<-(from+sum(mz1[[i]][[1]]$mass>=min_mz & mz1[[i]][[1]]$mass<=max_mz,na.rm=TRUE)-1)
				getpeaks[from:to,1]<-mz1[[i]][[1]]$mass[(mz1[[i]][[1]]$mass>=min_mz & mz1[[i]][[1]]$mass<=max_mz)]
				getpeaks[from:to,2]<-mz1[[i]][[1]]$intensity[(mz1[[i]][[1]]$mass>=min_mz & mz1[[i]][[1]]$mass<=max_mz)]				
				getpeaks[from:to,3]<-(mz1[[i]]$metaData$retentionTime)
				from<-(to+1)	  
			}
		}  
		if(progbar==TRUE){close(prog)};
		getpeaks[,4]<-seq(1,length(getpeaks[,1]),1)
		colnames(getpeaks)<-c("m/z","intensity","RT","measureID","partID","clustID","peakID")
		scans[[2]]<-getpeaks;
		rm(getpeaks);
		MSlist<-genMSlist();      
		MSlist[[4]]<-scans;
		rm(scans)
		MSlist[[1]][[1]]<-TRUE
		##########################################################################
		return(MSlist);
	      
}
