print.DocWordTable <-
function (x, file = NULL, sep = ";", ...) 
{
	res.DocWordTable <- x
	if (!inherits(res.DocWordTable, "DocWordTable")) stop("non convenient data")
    	cat("*The results are available in the following objects:\n\n")
        indice <-7
    	res <- array("", c(indice, 2), list(1:indice, c("name", "description")))
      res[1, ] <- c("$Ndoc", "Number of documents")
   	res[2, ] <- c("$Nlength", "Corpus size")
    	res[3, ] <- c("$Nword", "Vocabulary size")   	
    	res[4, ] <- c("$DocTerm", "Documents by Words table")
	res[5, ] <- c("$Tfreq", "Glossary")
    	res[6, ] <- c("$Nfreqword", "Frequencies Words")
	res[7, ] <- c("$Ndocword", "Frequencies of word in documents")       	 
    print(res[1:indice, ])
    	if (!is.null(file)) {
        	write.infile(res.DocWordTable, file = file, sep = sep)
        	print(paste("All the results are in the file", file))
    	}
}
