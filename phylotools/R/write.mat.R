#### Function add.mat as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010

write.mat <- 
function(supermat, file = NULL){
	if(is.null(file)){
	    stop("File name have to be specified.")
	}
    res <- dat2phy(supermat)	   
    writeLines(res, file)
}
