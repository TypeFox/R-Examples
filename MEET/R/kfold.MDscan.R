kfold.MDscan <-function(iicc,TF){

    require("seqinr")
	require("MEET")

    write.fasta <- get("write.fasta",pos="package:seqinr")
    read.fasta <- get("read.fasta",pos="package:seqinr")
	k<-length(iicc$DNAreal)
    	
	y<-read.fasta(file=TF,forceDNAtolower=FALSE)
	y<-y[-iicc$outsequence]
	write.fasta(sequences=y,names=c(1:length(y)),file.out="setTF.fa",open="w")
	validation.set_x <- iicc$DNAreal
   
    lenmotif<-iicc$lenmotif
    len_motif_inici<-(lenmotif-3)
    len_motif_final<-(lenmotif+3)
    num_motif<-iicc$nummotif
    direction<-iicc$direction
    call.MDscan<-iicc$MDscan
    
     Scoretotal <- lapply(seq(len_motif_inici,len_motif_final , 1), function(x){})
   
    for (m in c(1:length(y))) {
		 
		x<-y[-m]
		write.fasta(sequences=x,names=c(1:length(x)),file.out="factor.fa",open="w")
		training.set<-switch(iicc$alignment, "CLUSTALW"=align.clustalw(filein="factor.fa", fileout="background.fa", call=iicc$call.clustalw), "MUSCLE"=align.muscle(filein="factor.fa", fileout="background.fa", gapopen=iicc$gapopen, maxiters=iicc$maxiters, gapextend=iicc$gapextend, call=iicc$call.muscle),"MEME"=align.MEME(filein="factor.fa",fileout="background.fa",iicc),"NONE"=Read.aligned("factor.fa"), stop("Alignment method not included"))
		
		for (len_motif in seq(len_motif_inici,len_motif_final,1)){
	    
			resultats_MDscan<-run.read.MDscan(validation.set_x,k, len_motif, num_motif, call.MDscan)
			Score<-scoreMDscan(resultats_MDscan,k,training.set,direction)
			Scoretotal[[len_motif-len_motif_inici+1]]<-rbind( Scoretotal[[len_motif-len_motif_inici+1]],Score)  
	  }
	}
	SCOREtotal <- lapply(Scoretotal, function(x){as.vector(t(x))})
                
	return(SCOREtotal) 
	   
    }

