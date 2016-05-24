dataframe2fas <-
function(x, file = NULL)
{
    if(is.null(file)){
	"File name to save the results have to be given."
	}
	if(!is.data.frame(x)){
    x <- as.data.frame(x)
    }
    if(ncol(x)==1){
    x <- cbind(rownames(x), x)
    colnames(x) <- c("Seqnames", "squence")
    }
	if(! ncol(x) == 2){
	stop("Wrong dimention: the input dataframe must be \n either in 1 or 2 dimentions.")
	}
	dnanames <- as.character(x[,1])
	dnas <- as.character(x[,2])
	result <- c()
	for(i in 1:length(dnanames)) {
	   result1 <- dnanames[i]
	   result1 <- paste(">", result1, sep = "")
	   result2 <- dnas[i]
	   if(i == 1){
	    result <- c(result1, result2)
		}
	   if(i > 1){
        result <- c(result, result1, result2) 
       }
	}
	if(!is.null(file)){
	  writeLines(result, file)
	}
	result <- as.fasta(result)
	return(result)
}

