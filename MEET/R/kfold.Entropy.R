kfold.Entropy <-function(iicc,TF){

	require("MEET")
	require("seqinr")
	
	write.fasta <- get("write.fasta",pos="package:seqinr")
    read.fasta <- get("read.fasta",pos="package:seqinr")
    
    Prob<-as.numeric(iicc$background)
    iicc$missing.fun=TRUE
    k<- iicc$longDNA<-length(iicc$DNAreal)
	x<-read.fasta(file=TF,forceDNAtolower=FALSE)
    x<-x[-iicc$outsequence]
    write.fasta(sequences=x,names=c(1:length(x)),file.out="setTF.fa",open="w")
	
    	   print("factor")
    factor<-switch(iicc$alignment, "CLUSTALW"=align.clustalw(filein="setTF.fa", fileout="Sq.fa", call=iicc$clustalw), "MUSCLE"=align.muscle(filein="setTF.fa", fileout="Sq.fa", gapopen=iicc$gapopen, maxiters=iicc$maxiters, gapextend=iicc$gapextend, call=iicc$muscle),"MEME"=align.MEME(filein="setTF.fa",fileout="Sq.fa",iicc),"NONE"=Read.aligned("setTF.fa"), stop("Alignment method not included"))
	   print("factorfet")	
	validation.set_x <- iicc$DNAreal

    out<-lapply(seq(1, length(iicc$vector), 1), function(r){lapply(seq(1, length(x), 1), function(t){})})
       
	Resultats<-lapply(seq(1, length(iicc$vector), 1),function(i){
		
		q<-iicc$q<-iicc$correction_1rOrdre<-iicc$Herror<-iicc$HXmax<-iicc$Redundancia_corregida<-iicc$ErrorHX<-iicc$classentropy<-iicc$Entropy<-NULL
		
	    q<-iicc$q<-iicc$vector[i]
					  
		if (q==1) {iicc$classentropy<-"Shannon"
					  }else{
					iicc$classentropy<-"Renyi"  
					  }
					  
	    iicc$correction_1rOrdre <- correction.entropy(q,length(x),1,iicc)
	    iicc$Herror<-slot(iicc$correction_1rOrdre,"Herror")
		iicc$HXmax<-	 iicc$Herror[length(x)]
		iicc$Redundancia_corregida<-CalculRedundancy(factor,q,iicc)
	    iicc$ErrorHX<-slot(iicc$correction_1rOrdre,"sderror")[length(x)]
					  
		memory<-Hmemory(iicc,factor)
		iicc<-c(iicc,memory)

					  
		out[[i]]<-lapply(seq(1, length(x), 1), function(m){
	      
						 y<-x[-m]
						 write.fasta(sequences=y,names=c(1:length(x)),file.out="Sq.fa",open="w")
	  
						 training.set<-switch(iicc$alignment, "CLUSTALW"=align.clustalw(filein="Sq.fa", fileout="background.fa", call=iicc$call.clustalw), "MUSCLE"=align.muscle(filein="Sq.fa", fileout="background.fa", gapopen=iicc$gapopen, maxiters=iicc$maxiters, gapextend=iicc$gapextend, call=iicc$call.muscle),"MEME"=align.MEME(filein="Sq.fa",fileout="background.fa",iicc),"NONE"=Read.aligned("Sq.fa"), stop("Alignment method not included"))
	print("aliniat")
						 zout<-(sapply(X=c(1:(k-ncol(training.set)+1)),
						       FUN = function(X, training.set, validation.set_x,iicc) {
							   seq.rand <-validation.set_x[X:(X+ncol(training.set)-1)]
							   Hread(training.set.mes.rand=training.set,val.set=seq.rand,iicc)
							   }, training.set=training.set, iicc=iicc, validation.set_x=validation.set_x))
						 })
	  unlist(out[[i]])
	  
	})
	
	return(Resultats)
    }

