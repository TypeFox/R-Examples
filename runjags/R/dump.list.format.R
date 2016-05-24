dump.format <- function(namedlist=list(), checkvalid=TRUE, convertfactors=TRUE){
	
	data <- namedlist
	
	if(identical(data, list()))
		return('')
	
	if(length(data)==2 & is.null(names(data)) & class(data[[1]])=="character" & length(data[[1]])==1){  # allows old style dump.format (length = 1)
		names <- data
		data <- list(data[[2]])
		names(data) <- names[[1]]
	}
	
	if(inherits(data, 'data.frame'))
		data <- as.list(data)
	if(!inherits(data,"list") || length(data)==0) stop("Data must be provided as a named list or data frame", call.=FALSE)
	if(any(names(data)=="") || is.null(names(data))) stop("Data must be provided as a named list or data frame", call.=FALSE)
	if(length(unique(names(data)))!=length(data)) stop('All elements in the data list must have unique names', call.=FALSE)
	
	if(convertfactors){
		for(c in which(sapply(data,class)=='factor'))
			data[[c]] <- as.numeric(data[[c]])
	}
	
	if(checkvalid){
		valid <- checkvalidforjags(data)
		if(!valid$valid) stop(paste("The following problem was identified in the data provided:  ", valid$probstring, sep=""))				
	}
	
	variable = names(data)
	value <- data
	
	if(any(variable==".RNG.name")){
		n <- which(variable==".RNG.name")
		split <- strsplit(value[[n]], split="")[[1]]
		if(split[1]!="\"" & split[length(split)]!="\""){
			split <- c("\"", split, "\"")
			value[[n]] <- paste(split, collapse="")
		}
	}
	
	output.string <- ""
	for(i in 1:length(variable)){
		if(length(value[[i]])==1 && length(dim(value[[i]]))==0){
			value.string <- as.character(value[[i]])
		}else{
			dims <- dim(value[[i]])
			if(length(dims) > 1){
				value.string <- "structure(c("			
			}else{
				value.string <- "c("
			}
			value.string <- paste(value.string, paste(value[[i]], collapse=", "), ")", sep="")
			if(length(dims) > 1){
				value.string <- paste(value.string, ", .Dim = c(", paste(dims, collapse=", "), "))", sep="")
			}
		}
		output.string <- paste(output.string, "\"", variable[[i]], "\" <- ", value.string, "\n", sep="")
	}
		
	return(output.string)
}



list.format <- function(data=character(), checkvalid=TRUE){
	
	if(! inherits(data, c("character","runjagsdata","runjagsinits")) || length(data)==0) stop("Data must be provided as a character string in the R dump format")
	
	if(all(data=='' | data=='\n'))
		return(list())
	
	out <- vector('list', length=length(data))

	for(i in 1:length(data)){
		if(data[i]==""){
			out[[i]] <- list()
		}else{
			str <- data[i]
      # Remove anything following a comment:
      str <- gsub('#[^\n]*', '', str)
      # Remove leading and trailing white space:
			str <- gsub('^ *\n', '', str)
			str <- gsub('\n *\n$', '\n', str)
			
      str <- gsub("<-", "=", str)
			str <- gsub("`", "", str)
			str <- gsub("= \n", "=", str)
			str <- gsub("^\n", "", str)
			str <- gsub("\n\n", "\n", str)
			str <- gsub("\n\n", "\n", str)
			str <- gsub("\n\n", "\n", str)
      str <- gsub(",\n.Dim", ", .Dim", str, fixed=TRUE)
			str <- gsub("\n", ",", str)
			
			if(str!='' && strsplit(str, split="")[[1]][length(strsplit(str, split="")[[1]])] == ",") str <- paste(strsplit(str, split="")[[1]][1:(length(strsplit(str, split="")[[1]])-1)], collapse="")
			out[[i]] <- eval(parse(text=paste('list(', str, ')'))) 
		}
	}
	
	if(length(data)>1){
		names(out) <- paste('Chain.', 1:length(data), sep='')
	}else{
		out <- out[[1]]
	}
	
	if(checkvalid){
		valid <- checkvalidforjags(out)
		if(!valid$valid)
			stop(paste("The following problem was identified in the data provided:  ", valid$probstring, sep=""))				
	}
	
	return(out)
}
