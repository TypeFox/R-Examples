kfold.Divergence <-function(iicc,TF){

	require("MEET")
	require("seqinr")
	
    write.fasta <- get("write.fasta",pos="package:seqinr")
    read.fasta <- get("read.fasta",pos="package:seqinr")
    
    Prob<-iicc$background
    #iicc$missing.fun=TRUE
   	k<-iicc$longDNA<-length(iicc$DNAreal)
    x<-read.fasta(file=TF,forceDNAtolower=FALSE)
    x<-x[-iicc$outsequence]
    write.fasta(sequences=x,names=c(1:length(x)),file.out="setTF.fa",open="w")
	
	factor<-switch(iicc$alignment, "CLUSTALW"=align.clustalw(filein="setTF.fa", fileout="Sq.fa", call=iicc$clustalw), "MUSCLE"=align.muscle(filein="setTF.fa", fileout="Sq.fa", gapopen=iicc$gapopen, maxiters=iicc$maxiters, gapextend=iicc$gapextend, call=iicc$muscle),"MEME"=align.MEME(filein="setTF.fa",fileout="Sq.fa",iicc),"NONE"=Read.aligned("setTF.fa"), stop("Alignment method not included"))
	
    validation.set_x <- iicc$DNAreal
       
    out<-lapply(seq(1, length(iicc$vector), 1), function(r){lapply(seq(1, length(x), 1), function(t){})})
	
    
	Resultats<-lapply(seq(1, length(iicc$vector), 1),function(i){
					
		
		iicc$D <- iicc$correction_1rOrdre <- iicc$probparella<-NULL
		iicc$correction_1rOrdre_m1<-iicc$Mperfil<-iicc$q<- NULL
		iicc$HXmax<-iicc$HX<-iicc$q<-NULL
		iicc$Divergence<-iicc$interA<-iicc$interB<-iicc$classentropy<-NULL
			
		q<-iicc$q<-iicc$vector[[i]]
		
		if (q==1) {iicc$classentropy<-"Shannon"
					  }else{
				   iicc$classentropy<-"Renyi"  
					  }

		iicc <- detector_2nOrdre_init(factor, val.set=NULL, iicc)
		memory<-MImemory(iicc,factor)
		iicc<-c(iicc,memory)
			
 		out[[i]]<-lapply(seq(1, length(x), 1), function(m){	
						 
				y<-x[-m]
				write.fasta(sequences=y,names=c(1:length(x)),file.out="factor.fa",open="w")
						 
				training.set<-switch(iicc$alignment, "CLUSTALW"=align.clustalw(filein="factor.fa", fileout="background.fa", call=iicc$call.clustalw), "MUSCLE"=align.muscle(filein="factor.fa", fileout="background.fa", gapopen=iicc$gapopen, maxiters=iicc$maxiters, gapextend=iicc$gapextend, call=iicc$call.muscle),"MEME"=align.MEME(filein="factor.fa",fileout="background.fa",iicc),"NONE"=Read.aligned("factor.fa"), stop("Alignment method not included"))
										
				zout <- sapply( X=c(1:(iicc$longDNA-ncol(training.set)+1)),
						FUN = function(X, training.set, validation.set_x, iicc) {seq.rand <-validation.set_x[X:(X+ncol(training.set)-1)]
						MIread(training.set=training.set, val.set= seq.rand, iicc) 
						}, training.set=training.set, validation.set_x=validation.set_x, iicc=iicc)
			
	  })
	  unlist(out[[i]])
	})


 return(Resultats)
}

