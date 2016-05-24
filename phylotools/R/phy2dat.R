#### Function phy2dat as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010

phy2dat <-
function(x) {
	dd <- x
	firstl <- substring(dd[1], regexpr("[0-9]",dd[1]), nchar(dd[1]))
	nsp <- as.numeric(substring(firstl, 1,regexpr(" ", firstl)-1))
    nsite <- as.numeric(substring(firstl, regexpr(" ", firstl)+1, nchar(firstl)))
	x <- x[-1]
	if(nsp == length(x)){
	   col2 <- c()
	   col1 <- c()
	   for(i in 1:length(x)){
	      col1[i] <- substring(x[i], 1, regexpr(" ", x[i])-1)
		  col2[i] <- substring(x[i], regexpr(" ", x[i])+1, nchar(x[i]))
          col2[i] <- gsub(" ", "",col2[i])		  
	   }
	options(stringsAsFactors = FALSE)  
	col1 <- as.character(col1) 
	col2 <- as.character(col2)
	result <- data.frame(col1, col2)
	return(result)
	}
	else {
	    x <- x[grepl("[ATGC-]", x) > 0]
	    seqNam <- substr(x, 1, regexpr(" ", x) - 1)
	    seqNam <- seqNam[-which(seqNam == "")]
	    nspecies <- length(seqNam)
	    
	    nBlock <- length(x)/length(seqNam)
	    col.add <- rep(seqNam, nBlock)
	    
		## Core program adopted from Qiong Ding
	    for (i in 1:nspecies) {
	    	nBlock = length(x)/length(seqNam)
	    	rNam = ((1:nBlock) - 1) * nspecies + i
	    	stri = gsub("[ ]","", gsub(seqNam[i], "", paste(x[rNam],  collapse = "")))
	    	stri = toupper(stri)
	    		if (i == 1){ 
	    			DNA = stri
	    		}
	    		if (i > 1){ 
	    			DNA = c(DNA, stri)
	    		}	 
	    }
	    options(stringsAsFactors = FALSE)
	    DNA <- data.frame(seqNam,DNA)
        return(DNA)
    }	
}

