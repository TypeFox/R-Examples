##Recodes the variables of a data frame given a recoding specification
recode.variables<-function(data,recodes){
	recode.other<-function(var){
		if(is.factor(var)) stop("use recode.factor to recode factors")
		warning.flag<-TRUE
	    result <- var
		else.target<-""
		if(else.term!=""){
	        else.target <- eval(parse(text = strsplit(else.term, "->")[[1]][2]))
			result[1:length(var)] <- else.target
	    }
	    if(is.numeric(var)){
	        Lo <- min(var, na.rm = TRUE)
	        Hi <- max(var, na.rm = TRUE)
	    }else{
			Lo <-""
			Hi <-max(var, na.rm = TRUE)
		}
		for(term in recode.list){
			if(0 < length(grep(":", term))){
				if(is.character(var) && warning.flag){
					warning("Recoding a range of characters may not do what you think it does.\n Example: '15' is less than '9'.")
					warning.flag<-FALSE
				}
				range <- strsplit(strsplit(term, "->")[[1]][1], ":")
				low <- eval(parse(text = range[[1]][1]))
				high <- eval(parse(text = range[[1]][2]))
				if(high<low) next 
				target <- eval(parse(text = strsplit(term, "->")[[1]][2]))
				result[(var >= low) & (var <= high)] <- target
			}else{
				set <- eval(parse(text = strsplit(term, "->")[[1]][1]))
				target <- eval(parse(text = strsplit(term, "->")[[1]][2]))
		        for (val in set) {
					if (is.na(val)) 
				 		result[is.na(var)] <- target
		            else{ 
						result[var == val] <- target
					}	
				}
			}
		}
		return(result)
	}

	recode.factor<-function(var){
		if(!is.factor(var)) stop("var must be a factor")
		result<-var
		else.target<-""
		if(else.term!=""){
	        else.target <- eval(parse(text = strsplit(else.term, "->")[[1]][2]))
			if(!(else.target %in% levels(result))){
				levels(result)<-c(levels(result),else.target)
			}
			result<-factor(rep(else.target,length(var)),levels=else.target)
	    }

		for(term in recode.list){
			Lo<-levels(var)[1]
			Hi<-levels(var)[length(levels(var))]
			if(0 < length(grep(":", term))){
				range <- strsplit(strsplit(term, "->")[[1]][1], ":")
				low <- eval(parse(text = range[[1]][1]))
				low<-which(levels(var)==low)[1]
				if(is.na(low)) stop(paste("Lower value in range not a valid factor level.",term))
				high <- eval(parse(text = range[[1]][2]))
				high <- which(levels(var)==high)[1]
				if(is.na(high)) stop(paste("upper value in range not a valid factor level.",term))
				if(high<low) stop(paste("Upper value must be ordered after lower value in the factor ordering.",term))
	
				target <- eval(parse(text = strsplit(term, "->")[[1]][2]))
				set<-levels(var)[low:high]
				if(!(target %in% levels(result))){
					levels(result)<-c(levels(result),target)	
				}
				result[var %in% set] <- target
				set<-setdiff(set,target)
				levels(result)<-ifelse(levels(result) %in% set,NA,levels(result))
			}else{
				set <- eval(parse(text = strsplit(term, "->")[[1]][1]))
				target <- eval(parse(text = strsplit(term, "->")[[1]][2]))
		        for (val in set) {
					if(!(target %in% levels(result))){
						levels(result)<-c(levels(result),target)
					}
					if (is.na(val)) 
				 		result[is.na(var)] <- target
		            else{ 
						result[var == val] <- target
						if (!is.na(val) && !is.na(target) && val != target){
							levels(result)<-ifelse(levels(result)==val,NA,levels(result))
						}
					}	
				}
			}
		}
		return(result)
	}

	if(!is.data.frame(data)) data<-as.data.frame(data)
	recode.list <- strsplit(recodes, ";")[[1]]
	else.term<-""
	else.ind<-c()
	for(i in 1:length(recode.list)){
		first.part<-strsplit(recode.list[[i]],"->")[[1]][1]
		if(length(grep("else",first.part))>0 && length(grep("'",first.part))<1){
			else.term<-recode.list[[i]]
			else.ind<-c(else.ind,-i)
		}
	}
	if(length(else.ind)>0) recode.list<-recode.list[else.ind]	
	result.data<-data.frame(1:dim(data)[1])
	for(variable in data){
		if(is.factor(variable)){
			result.data<-data.frame(result.data,recode.factor(variable),stringsAsFactors=FALSE)
		}else result.data<-data.frame(result.data,recode.other(variable),stringsAsFactors=FALSE)
	}
	return(result.data[-1])
} 