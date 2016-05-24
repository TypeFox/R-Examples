#### Function seq2fasta as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010

seq2fasta <- 
function(file){
    fil <- suppressWarnings(readLines(file))
    seq <- fil[max(which(grepl("\\^\\^", fil)))+1]
    nam <- substr(file, start = max(unlist(gregexpr("/", file)))+1, stop = nchar(file))
    nam1 <- paste(">",nam, sep = "")
    res <- c(nam1, seq )
    newfilnam <- gsub("\\.seq", "\\.fasta", nam)
    writeLines(res, newfilnam)
	return(res)
}
