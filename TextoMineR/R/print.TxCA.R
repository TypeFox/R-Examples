print.TxCA <-
function (x, file = NULL, sep = ";", ...) 
{
	res.TxCA <- x
	if (!inherits(res.TxCA, "TxCA")) stop("non convenient data")
	cat("**Results for AC and Aggregated Lexical Table (TxCA)**\n")
    	cat("*The results are available in the following objects:\n\n")
       indice <-11 
    	res <- array("", c(indice, 2), list(1:indice, c("name", "description")))
   	res[1, ] <- c("$TableSummary", "Summary of the information")
    	res[2, ] <- c("$Glossary", "Glossary of words") 
        res[3, ] <- c("$Table", "Table was analyzed by CA")  
    	res[4, ] <- c("$DocTermR", "Table Documents by Words")
    	res[5, ] <- c("$res.agg", "Result of aggregation")
    	res[6, ] <- c("$Nfreqword", "Frequencies of words")
	   res[7, ] <- c("$Ndocword", "Frequencies of words in documents")
    	res[8, ] <- c("$res.ca", "Results of correspondence analysis")
	res[9,] <- c("$VCr","Cramers V")
	res[10,] <- c("$Inertia","Total inertia")
	res[11,] <- c("$res.meta", "Results of Metakeys and Metadocs")
    print(res[1:indice, ])
    	if (!is.null(file)) {
        	write.infile(res.TxCA, file = file, sep = sep)
        	print(paste("All the results are in the file", file))
    	}
}
