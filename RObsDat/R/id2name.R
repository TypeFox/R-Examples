# convert a numeric to string

id2name <- function(dataframe){
	dataframe.new <- dataframe
	names.alt <- sub("id","",names(dataframe))
	for(i in grep("id",names(dataframe))){
		if(names.alt[i]=="version") next #no conversion for version ID
		repl <- getMetadata(names.alt[i], ID=dataframe[,i], EXACT=TRUE)

		if(!is.null(repl$Name)){
		     dataframe.new[,i] <- repl$Name
		}  else if (!is.null(repl$Description)){
		     dataframe.new[,i] <- repl$Description
		} else {
		     dataframe.new[,i] <- repl$ID
		}
	}
	if(NCOL(dataframe.new)!=NCOL(dataframe)){
		stop("lost columns in id2name")
	}
	names(dataframe.new) <- names.alt
	return(dataframe.new)
}
