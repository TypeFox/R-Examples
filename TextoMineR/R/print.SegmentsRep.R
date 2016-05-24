print.SegmentsRep <-function (x, file = NULL, sep = ";", ...) 
{
	res.segments <- x
	if (!inherits(res.segments, "SegmentsRep")) stop("non convenient data")
	cat("**Results for Aggregated repeated segments (segments)**\n")
    	cat("*The results are available in the following objects:\n\n")
    	indice <-4 
    	res <- array("", c(indice, 2), list(1:indice, c("name", "description")))
   	res[1, ] <- c("$nbseg", "Number of repeated segments")
    	res[2, ] <- c("$tab.seg", "Table Documents by segments")
    	res[3, ] <- c("$list.tot.segment", "complete segments's list")
    	res[4, ] <- c("$Glossary.segments", "Glossary of repeated segments")
    print(res[1:indice, ])
    	if (!is.null(file)) {
        	write.infile(res.segments, file = file, sep = sep)
        	print(paste("All the results are in the file", file))
    	}
}
