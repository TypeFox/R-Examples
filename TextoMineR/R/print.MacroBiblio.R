print.MacroBiblio <-
function (x, file = NULL, sep = ";", ...) 
{
	res.MacroBiblio <- x
	if (!inherits(res.MacroBiblio, "MacroBiblio")) stop("non convenient data")
	cat("**Results of Analysis of bibliography  (MacroBiblio)**\n")
    	cat("*The results are available in the following objects:\n\n")
	indice <-12 
    	res <- array("", c(indice, 2), list(1:indice, c("name", "description")))
   	res[1, ] <- c("$Corpus", "Summary of the information about corpus")
    	res[2, ] <- c("$Glossary", "Glossary of words") 
    	res[3, ] <- c("$DocTermR", "Documents by words table")
        res[4, ] <- c("$Tagreg", "Lexical aggregate table")
        res[5, ] <- c("$Metakeys.Metadocs", "Representation of words/documents with higher contribution ")
    	res[6, ] <- c("$res.CA", "Results of  correspondence analysis direct")
    	res[7, ] <- c("$res.CA.Agreg", "Result of Correspondence analysis aggregate")
	res[8, ] <- c("$CharWord", "characteristic words in each group of the aggregation variable")
    	res[9, ] <- c("$res.CHCPC", "Result of constrained hierarchical clustering")
	res[10,] <- c("$res.MFACT","Result of multiple factor analysis for contingency tables")
	res[11,] <- c("$OrdWord", "words and their coordinates in the first dimension")
        res[12,] <- c("$pioneers", "pioneers articles")     
    print(res[1:indice, ])
    	if (!is.null(file)) {
        	write.infile(res.MacroBiblio, file = file, sep = sep)
        	print(paste("All the results are in the file", file))
    	}
}
