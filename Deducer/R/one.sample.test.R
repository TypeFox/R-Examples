one.sample.test<-function(variables,data=NULL,test=t.test,...){
	arguments <- as.list(match.call()[-1])
	vars<-eval(substitute(variables),data, parent.frame())
	if(length(dim(vars))<1.5){
		vars<-d(vars)
		fn<-arguments$variables
		names(vars)<-if(is.call(fn)) format(fn) else as.character(fn)
	}
	data <- vars
	tests<-list(NULL)
	for(i in 1:ncol(data)){
		if(is.character(data[[i]])){
			cat(paste("'",names(data)[i],
					"' is a character vector. Attempting to coerse into a numeric one.\n",sep=""))
			data[[i]]<-as.numeric(data[[i]])
		}
		if(is.factor(data[[i]])){
			cat(paste("'",names(data)[i],
					"' is a factor. Attempting to coerse into numeric.\n",sep=""))
			data[[i]]
			tmp<-as.numeric(as.character(data[[i]]))
			if(length(na.omit(tmp))==length(na.omit(data[[i]])))
				data[[i]] <- tmp
			else
				data[[i]] <- as.numeric(data[[i]])
		}
		if(is.logical(data[[i]])){
			data[[i]]<-as.numeric(data[[i]])
		}
		if(!is.numeric(data[[i]])){
			warning(paste("'",names(data)[i],
					"' is not numeric. It will be dropped.\n",sep=""))
			next
		}
		if(length(na.omit(data[[i]]))<3){
			warning(paste("'",names(data)[i],
					"' has fewer than 3 valid values. It will be dropped.",sep=""))
			next
		}
		tests[[names(data)[i]]]<-test(na.omit(data[[i]]),...)
	}
	multi.test(tests)
}