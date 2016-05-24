#### Delete the duplicated sequence in fasta file, 
#### if duplicated sequences have been found, only the first one will be retained.
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Oct- 30- 2011
#### Bug fixed  Nov.-01-2010

uniquefasta <- 
function(fasta){
    names.fas <- gnames.fas(fasta)
    times <- as.data.frame(table(names.fas))
    duplctd <- as.character(times[,1])[times[,2] > 1]
	if(!any(times[,2] > 1)){
	    cat("no duplicate names found.\n")
		return(fasta)
	}
	res0 <- fasta
    for(i in 1:length(duplctd)){
        index.name <- which(duplctd[i] == gnames.fas(as.fasta(res0)))
        index.seq <- index.name*2
        index <- sort(c(index.seq -1, index.seq))
		index.del <- index[-c(1,2)]
        res0 <- res0[-index.del]
    }
    res <- as.fasta(res0)
	return(res)
}
