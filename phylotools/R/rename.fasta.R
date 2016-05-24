#### Renaming of the sequences in a fasta file, given a reference table
#### This function is part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Oct- 30- 2010


rename.fasta <- 
function(fas, ref, fil = NULL, prefix = NULL){
    ## Convert the input ref to dataframe
    if(!is.data.frame(ref)){
 	    ref <- as.data.frame(ref)    
	}
	
	### if(!length(dim(ref)) == 2){
	###     stop("The input ref table must be a dataframe")
	### }
	
    names.fas <- gnames.fas(fas)
    ref[,1] <- paste(prefix, ref[,1], sep = "")
    names.ref <- ref[,1]
    names.sub <- as.character(ref[,2])
    new.names.fas <- rep(NA, length(names.fas))
    
	##Which name to apply
    for(i in 1:nrow(ref)){
        for(j in 1:length(names.fas)){
    	    if(names.fas[j] == names.ref[i]){
    		    new.names.fas[j] <- names.sub[i]
    		}
    	}
    }
	
	##Give out Warning if any of the names in fasta can not be found in the ref table
    if(any(is.na(new.names.fas))){
        unfind <- names.fas[which(is.na(new.names.fas))]
        warning(paste(paste(unfind, collapse = ", "),"can not be found in the refence table"))
    }
	
	## Rename
    res <- rename.fas(fas, paste(">",new.names.fas, sep = ""))
    
	## If specified the file name, write to the file.
	if(!is.null(fil)){
	    write.fasta(res, fil)
	}
	return(res)
}

