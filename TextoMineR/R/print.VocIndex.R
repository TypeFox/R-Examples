print.VocIndex <-
function (x, file = NULL, sep = ";", ...) 
{
	res.VocIndex<- x
	if (!inherits(res.VocIndex, "VocIndex")) stop("non convenient data")
	cat("**Induce of Vocabulary  (VocIndex)**\n")
      	indice <-3 
    	cat("*The results are available in the following objects:\n\n")
    	res <- array("", c(indice, 2), list(1:indice, c("name", "description")))
   	res[1, ] <- c("$RegVoc", "Regular or stable vocabulary")
    	res[2, ] <- c("$LocalVoc","Local or specialized vocabulary")
    	res[3, ] <- c("$VocIndex","Vocabulary index table  ")	
    print(res[1:indice, ])
    	if (!is.null(file)) {
        	write.infile(res.VocIndex, file = file, sep = sep)
        	print(paste("All the results are in the file", file))
    	}
}
