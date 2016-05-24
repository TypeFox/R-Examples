#### Split the fasta file to several files, according to the group each sequence belonging to.
#### This function is part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Oct- 30- 2011


fasta.split <- 
function(fasta, ref, save2disk = FALSE){
    
    ## if(!is.integer(length(fasta)/2)){
    ##    stop("make sure the input data contains sequences in fasta format.")
    ## }
    
    fas.nams <- gsub(">", "", fasta[seq(1, length(fasta), by = 2)])
    
    names.ref <- as.character(ref[,1])
    level.names <- as.character(ref[,2])
    level.names.unique <- unique(level.names)
    res0 <- list()
    fasta.group <- rep(NA, length(fasta)/2)
	## Find the group number for each sequence
    for(i in 1:(length(fasta)/2)){
        for(j in 1:(nrow(ref))){
    	    if(fas.nams[i] == names.ref[j]){
    		    fasta.group[i] <- level.names[j]
    		}
    	}
    }
    
	## Save the sequence to the list and to the disk
    for(k in 1:length(level.names.unique)){
        seqnams.ind <- seq(1, length(fasta), by = 2)[fasta.group == level.names.unique[k]]
    	res.nam <- gsub(">", "",fasta[seqnams.ind])
    	seq.ind <- seq(1, length(fasta), by = 2)[(fasta.group == level.names.unique[k])] + 1
    	res.seq <- fasta[seq.ind]
    	res0[[k]] <- dataframe2fas(data.frame(res.nam, res.seq))
    	if(save2disk){
		    write.fasta(res0[[k]], paste("group",k,".fasta", sep = ""))
		}
    }
    return(res0)
}

