
seqsep <- function (seqdata,sl=1, sep="-") {
	for (i in 1:length(seqdata)) {
		oseq <- seqdata[i]
		seql <- nchar(oseq)

		if ((seql %% sl)==0) {
			nbs <- seql/sl
			if (nbs>0) subseq <- substr(oseq,1,sl)
			if (nbs>1) {
				for (j in 2:nbs) {
					start <- ((j-1)*sl)+1
					stop <- start+(sl-1)							
					subseq <- paste(subseq,sep,substr(oseq,start,stop),sep="")
					}
				}
					seqdata[i] <- subseq
				}
			else 
				stop("Number of characters does not match number of states*states length in sequence",i)
		}	
	return(seqdata)
}
