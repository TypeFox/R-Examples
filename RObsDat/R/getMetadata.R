getMetadata <- function(table, EXACT=FALSE, ...){
	#check for valid table name
	table.list <- c(CVtables(), "SpatialReference", "Site", "Method","Qualifier","QualityControlLevel","Sample","Source","Variable","OffsetType","Units","ISOMetadata", "LabMethods")
	if(!table %in% table.list){
		lower.match <- tolower(table.list) %in% table
		if(any(lower.match)){
			table <- table.list[lower.match]
		} else {
			stop("Undefined table ", table, " Valid values are: ", paste(table.list, collapse=", "))
		}
	}

	#query data
	entry <- NULL
	extras <- list(...)
	length.extras <- sapply(extras, length)
	if(length(extras) > 1 &  any(length.extras != 1)){
		todo("treating multiple additional arguments to getMetadata")
		browser()
	}
	if(length(extras) == 1 & EXACT){
		u.entries <- unique(extras[[1]])
		for(i in seq(along=u.entries)){
		    command <- paste('entry <- Iget',table,'(options("odm.handler")[[1]], ',names(extras),'="',u.entries[i],'", exact=TRUE )', sep='')
		    eval(parse(text=command))
		    if(!exists("meta")){
			    meta <- as.data.frame(matrix(NA, nrow=length(extras[[i]]), ncol=NCOL(entry)))
			    names(meta) <- names(entry)
		    }
		    if(NROW(entry)==0) stop(paste("No meta data found in table ", table," for ", names(extras), "==", u.entries[i]))
		    if(NROW(entry)>1) stop(paste("Multiple datasets found in table ", table," for ", names(extras), "==", u.entries[i], "Searching for a variable not part of the table?"))
		    #ToDo: Warning generated when accessing VariableNames and ID
		    meta[extras[[1]]==u.entries[i],] <- entry
		}
		return(meta)

	} else {
		command <- paste('entry <- Iget',table,'(options("odm.handler")[[1]], ...)', sep='')
		eval(parse(text=command))
		return(entry)
	}

}
