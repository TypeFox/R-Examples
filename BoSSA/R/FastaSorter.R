`FastaSorter` <-
function(Y)
{
	alignbrut <- read.dna(Y,format="fasta")
	align <- rename(alignbrut)
	InfoFasta(alignbrut,Y)
	write.dna(align,paste(substring(as.character(Y),1,5),"Sort.fas",sep=""),format="fasta",colsep="")
	align
}

