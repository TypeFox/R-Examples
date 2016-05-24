PredictMEME <-function(iicc){
    
    require("seqinr")
   
    read.fasta <- get("read.fasta",pos="package:seqinr")
	
	call.mast<-iicc$call.mast
    k<-length(iicc$DNA[[1]])
    num_motif<-iicc$nummotif
    len_motif<-iicc$lenmotif
    threshold<-iicc$threshold
    output<-lapply(seq(1,length(iicc$DNA), 1), function(x){})
    
      for(i in c(1:length(iicc$DNA))){
      
	    Sequence <- iicc$DNA[[i]]
	    write.fasta <- get("write.fasta",pos="package:seqinr")
	    write.fasta(Sequence, names="sequenciaEstudi", nbchar = k, file.out="sequenciaEstudi.fa",open="w")
	    system(paste(paste(paste(call.mast, "memeout/meme.txt -d sequenciaEstudi.fa -text -mt",sep=" "), threshold, sep=" "), "-hit_list", sep=" "))
	    meme_thresholds<-read.mast(paste("mast.meme.txt.sequenciaEstudi.fa.mt",threshold,sep=""),factor,k )
	    meme_thresholds[is.na(meme_thresholds)]<-1
	    meme_thresholds<-(as.numeric(meme_thresholds))
	  
	    DetectedFactors<-cbind(which(meme_thresholds!=1),meme_thresholds[meme_thresholds!=1])
	    if(nrow(DetectedFactors)==0){
	      output<-"No Binding Sites Found"
	    }else{
	      output[[i]]<-lapply(c(1:nrow(DetectedFactors)),function(x){cbind(Sequence=paste(Sequence[DetectedFactors[x,1]:(DetectedFactors[x]+len_motif)],sep="",collapse=""),pvalue=DetectedFactors[x,2], position=DetectedFactors[x,1])})
	    }
	    return(output)
    }

}

