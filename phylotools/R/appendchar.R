#### Function appendchar as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010

appendchar <-
function(mat, pattern = "?"){
    if(!(nchar(pattern) == 1)){
	    stop("only one character is allowed.")
	}
    if(any(!nchar(mat[,2]) == max(nchar(mat[,2])))){
        for(i in 1:length(mat[,2])){
            if(!nchar(mat[i,2])== max(nchar(mat[,2]))){ 
               mat[i, 2] <- paste(mat[i, 2], 
                          paste(rep( x = pattern, (max(nchar(mat[,2])) 
						  - nchar(mat[i,2]))), collapse = ""), sep = "")
            }
	    }
	}
	return(mat)
}

