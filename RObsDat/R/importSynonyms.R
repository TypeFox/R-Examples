importSynonyms <- function(file){
	datf <- read.table(file, header=TRUE, as.is=TRUE)
	if(paste(names(datf), collapse="") != "phrasetablekey"){
		stop("file seems not to be a synonym export")
	}
	for(i in 1:NROW(datf)){
		cat("Importing: ", datf$phrase[i], "\n")
		try(addSynonym(table=datf$table[i], phrase=datf$phrase[i], id=datf$key[i]))
	}
	cat(NROW(datf), " synonyms imported \n")
}

