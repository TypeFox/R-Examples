align.muscle <-function(filein, fileout="Sq.fa", gapopen=gapopen, maxiters=maxiters, gapextend=gapextend, call){
    read.fasta <- get("read.fasta", pos = "package:seqinr")
	strtocall <- paste(call, ' -in ',filein ,' -out',fileout,
				'-maxiters',maxiters,'-diags1',
				'-gapopen',gapopen)

	
	system(strtocall)
	sq<-read.fasta(file=fileout)
	sqalineada<-matrix(,length(sq),length(sq[[1]]))
	for (i in c(1:length(sq))){
		sqalineada[i,]<-sq[[i]]
	}
	return(sqalineada)
}

