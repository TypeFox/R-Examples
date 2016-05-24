print.MacroCaHcpc <-
function (x, file = NULL, sep = ";", ...) 
{
    res<- x
    if (!inherits(res, "MacroCaHcpc")) 
        stop("non convenient data")
    cat("**Results for the MacroCaHcpc **\n")
    indice <- 5
    res <- array("", c(indice, 2), list(1:indice, c("name", "description")))
	res[1, ] <- c("$Corpus", "Description of corpus")
	res[2, ] <- c("$res.TxCA", "Results of correspondence analysis")
    res[3, ] <- c("$res.TxCharClust", "characteristic documents and words of the clusters")
    res[4, ] <- c("$res.hcpc", "Results of hierarchical clustering")	
    res[5, ] <- c("$ncp", "number of dimensions preserved")
    print(res[1:indice, ])
    if (!is.null(file)) {
        write.infile(MacroCaHcpc, file = file, sep = sep)
        print(paste("All the results are in the file", file))
    }
}




    
      