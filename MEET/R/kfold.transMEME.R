kfold.transMEME <-function(iicc,TF) {
 
   require("MEET")
   require("seqinr")
    
    
    write.fasta <- get("write.fasta",pos="package:seqinr")
    read.fasta <- get("read.fasta",pos="package:seqinr")
      
   	k<-length(iicc$DNAreal)
	
	call.transfac2meme<-iicc$transfac2meme
	call.mast<-iicc$mast
	x<-read.fasta(file=TF)
	x<-x[-iicc$outsequence]
	write.fasta(sequences=x,names=c(1:length(x)),file.out="setTF.fa",open="w")
	factor<-switch(iicc$alignment, "CLUSTALW"=align.clustalw(filein="setTF.fa", fileout="Sq.fa"), "MUSCLE"=align.muscle(filein="setTF.fa", fileout="Sq.fa", gapopen=iicc$gapopen, maxiters=iicc$maxiters, gapextend=iicc$gapextend),"MEME"=align.MEME(filein="setTF.fa",fileout="Sq.fa",iicc),"NONE"=Read.aligned(TF), stop("Alignment method not included"))
	
    validation.set_x <- iicc$DNAreal
    
    meme_thresholds<-motiu<- vector("numeric",length=((k-ncol(factor)+1)*nrow(factor)))
    
    write.fasta(validation.set_x, names="sequenciaEstudi", nbchar = k, file.out="sequenciaEstudi.fa",open="w")
    
    threshold<-0.1
    
    for (m in c(1:nrow(factor))) {
		    
	    training.set<-factor[-m,]
	    
	    writeMEME(training.set,m)
	    
	    system(paste("./transfac2meme -bg bfile -pseudo 0.1 motif.dat>motif.meme"))
	    
	    system(paste(paste("./mast  motif.meme -d sequenciaEstudi.fa -text -mt", threshold, sep=" "), "-hit_list", sep=" "))
	    meme_thresholds[(((k-ncol(factor)+1)*(m-1)+1):((k-ncol(factor)+1)*m))]<-read.mast(paste("mast.motif.meme.sequenciaEstudi.fa.mt",threshold,sep=""),factor,k )
	    meme_thresholds<-(as.numeric(meme_thresholds))
	    motiu[(((k-ncol(factor)+1)*(m-1)+1):((k-ncol(factor)+1)*m))]<-motif.mast(paste("mast.motif.meme.sequenciaEstudi.fa.mt",threshold,sep=""),factor,k,m)
	    motiu<-(as.numeric(motiu))
	    
	}
	    
	
	w <- vector(mode='numeric',length=(k-ncol(factor)+1)*nrow(factor))
	for (i in c(1: nrow(factor))){
	      w[(iicc$position-round(ncol(factor)/2))-1+(k-ncol(factor)+1)*(i-1)+order(meme_thresholds[iicc$position+ c(-round(ncol(factor)/2):round(ncol(factor)/2))+(k-ncol(factor)+1)*(i-1)])[1]]<-1 
	    }
      
	MEME_logthresholds<-(-log10(as.numeric(meme_thresholds)))
	print(summary(MEME_logthresholds))
  	threshold<-sort(meme_thresholds[w==1])
  	s<-length(threshold)
  	threshold_final<-sapply(c(1:s),function(i){	
  		if (i==1){0.5*threshold[i]
  		}else{
   		(threshold[i]+threshold[i-1])*0.5}
  		})
  	threshold_final[s+1]<-(threshold[s]*2)
#   
  	MEME_logthresholds<-MEME_thresholds<- lapply(c(1:(s+1)), function(x){rep(0,nrow(factor)*(k-ncol(factor)+1))})	
  	area_MEME<-vector("numeric",length=length(MEME_thresholds))
  	perf<-lapply(seq(1, length(MEME_thresholds), 1), function(x){}) 
  	
 	    
     for (m in c(1:nrow(factor))) {
 		    
 	    training.set<-factor[-m,]
 	    writeMEME(training.set,m)
 	    write.fasta(validation.set_x, names="sequenciaEstudi", nbchar = k, file.out="sequenciaEstudi.fa",open="w")
 	    system(paste("./transfac2meme -bg bfile -pseudo 0.1 motif.dat>motif.meme"))
 	    
 	     for (i in c(1:(s+1))) {
  		      system(paste(paste("./mast  motif.meme -d sequenciaEstudi.fa -text -mt",threshold_final[i], sep=" "), "-hit_list", sep=" "))
  		      MEME_thresholds[[i]][(((k-ncol(factor)+1)*(m-1)+1):((k-ncol(factor)+1)*m))]<-read.mast(paste("mast.motif.meme.sequenciaEstudi.fa.mt",threshold_final[i],sep=""),factor,k )
  		      MEME_logthresholds[[i]]<-(-log10(as.numeric(MEME_thresholds[[i]])))
  		      }
  	}
 	

	return(MEME_logthresholds)
 }

