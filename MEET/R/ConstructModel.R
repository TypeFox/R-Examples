ConstructModel<-function(iicc,TF){
	require("MEET")
	require("seqinr")
	#require("fields")

    write.fasta <- get("write.fasta",pos="package:seqinr")
    read.fasta <- get("read.fasta",pos="package:seqinr")
	
	x<-read.fasta(file=TF,set.attributes=FALSE, forceDNAtolower=FALSE)
	iicc$length<-length(x)

    Rdetector<-lapply(seq(1, iicc$length, 1), function(sequence){
					  
					  iicc$outsequence<-sequence
	     
					  iicc$DNAreal<-c(iicc$DNA[[1]][1:(iicc$position-1)],x[[iicc$outsequence]], iicc$DNA[[1]][(iicc$position+length(x[[sequence]])):length(iicc$DNA[[1]])])
    
					  Rdec<-switch(iicc$method, "Entropy"=kfold.Entropy(iicc,TF=TF), "Divergence"=kfold.Divergence(iicc,TF=TF),"Qresiduals"=kfold.PCA(iicc,TF=TF) ,"MEME"=kfold.MEME(iicc,TF=TF),"MDscan"= kfold.MDscan(iicc,TF=TF),"MATCH"=kfold.MATCH(iicc,Seqin=TF),"transMEME"=kfold.transMEME(iicc,TF=TF),stop("Detection method not included"))
					 
					  Rdec
			 })
        
	p<-standardout(Rdetector,iicc)
	
	Scores<-p$Scores; labels<-p$labels
	
	ROCcurve<-ROCmodel(Scores, labels, iicc)
	
	return(ROCcurve)
	
}
