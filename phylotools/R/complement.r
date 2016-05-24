#### Function complement as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010


complement <- 
function(fas){
    if(!class(fas) == "fasta"){
	    stop("the input data must be in fasta format")
	}
    seqind <- seq(2, length(fas), by = 2)
	for(i in seqind){
	    fas[i] <- revComp(fas[i])
	}
    return(fas)
}

