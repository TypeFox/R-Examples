kfold.MEME <-function(iicc,TF){

    require("seqinr")
        
    write.fasta <- get("write.fasta",pos="package:seqinr")
    read.fasta <- get("read.fasta",pos="package:seqinr")
     
	call.meme<-iicc$meme
	call.mast<-iicc$mast
    k<-length(iicc$DNAreal )
    x<-read.fasta(file=TF,forceDNAtolower=FALSE)
	x<-x[-iicc$outsequence]
	    
    validation.set_x <- iicc$DNAreal
    num_motif<-iicc$nummotif
    len_motif<-iicc$lenmotif
    threshold<-0.1
    
    meme_thresholds<-lapply(seq(1,length(x) , 1), function(x){})
    write.fasta(validation.set_x, names="sequenciaEstudi", nbchar = k, file.out="sequenciaEstudi.fa",open="w")
    w<-NULL
  
    for (m in c(1:length(x))) {
    
		y<-x[-m]
		write.fasta(sequences=y,names=c(1:length(x)),file.out="background.fa",open="w")
		
		system(paste(paste(paste(paste(paste(call.meme,TF,sep=" "),"-oc memeout -dna -mod anr -w",sep=" "), len_motif,sep=" "), "-nmotifs", sep=" "),num_motif,sep=" "))
	   
	    system(paste(paste(paste(call.mast,"memeout/meme.txt -d sequenciaEstudi.fa -text -mt",sep=" "), threshold, sep=" "), "-hit_list", sep=" "))
	    meme_thresholds[[m]]<-read.mast(paste("mast.meme.txt.sequenciaEstudi.fa.mt",threshold,sep=" "),iicc$Transcriptionfactor,k)
	    meme_thresholds[[m]]<-(as.numeric(meme_thresholds[[m]]))
	    meme_thresholds[[m]][which(is.na(meme_thresholds[[1]]))]<-1
		ww<-vector(mode="numeric",length=length(meme_thresholds[[m]]))
	    
	    for (p in c(1:length(iicc$position))) {
	      for (j in c(-(ncol(iicc$Transcriptionfactor)):(ncol(iicc$Transcriptionfactor)))) {
		  if (meme_thresholds[[m]][iicc$position[p]+ j]!=1) {
		    ww[iicc$position[p]+ j]<-1}
		  }
	    }
	    w<-rbind(w,ww)
	}
	
	meme_thresholds<-unlist(meme_thresholds)
	w<-as.vector(t(w))
	
	    
 	threshold<-sort(meme_thresholds[w==1])
 
 	s<-length(threshold)
 	threshold_final<-sapply(c(1:s),function(i){	
 		if (i==1){0.5*threshold[i]
 		}else{
  		(threshold[i]+threshold[i-1])*0.5}
 		})
 	threshold_final[s+1]<-(threshold[s]*2)
 	
 	MEME_logthresholds<-MEME_thresholds<- lapply(c(1:(s+1)), function(g){ lapply(seq(1,length(y) , 1), function(h){})})	
 		
	MEME_logthresholds<-lapply(c(1:(s+1)), function(i){

			MEME_thresholds[[i]]<-lapply(seq(1,length(x) , 1), function(m){
		  
			  y<-x[-m]
			  write.fasta(sequences=y,names=c(1:length(x)),file.out="background.fa",open="w")
 		      system(paste(paste("./mast  memeout/meme.txt -d sequenciaEstudi.fa -text -mt",threshold_final[i], sep=" "), "-hit_list", sep=" "))
			  zout<-read.mast(paste("mast.meme.txt.sequenciaEstudi.fa.mt",threshold_final[i],sep=""),iicc$Transcriptionfactor,k )
			  return(zout)
		      })
							
	out<-(-log10(as.numeric(unlist(MEME_thresholds[[i]]))))
	return(out)
	
 	})

	return(MEME_logthresholds)
	
  }

