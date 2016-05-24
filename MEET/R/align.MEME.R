align.MEME <-function(filein,fileout="Sq.fa",iicc){
	require("MEET")
	num_motif<-iicc$nummotif
	len_motif<-iicc$lenmotif
	call.meme=iicc$call.meme
	system(paste(paste(paste(paste(call.meme,filein,sep=" "),"-oc memeout -dna -mod anr -w", len_motif,sep=" "), "-nmotifs", sep=" "),num_motif,sep=" "))
	
	resultat<-readLines(paste(getwd(),"/memeout/meme.txt", sep="", collapse=NULL))
	y <- readMEME(resultat, num_motif)

	Seqalineada<-y[[1]]
	return(Seqalineada)
}

