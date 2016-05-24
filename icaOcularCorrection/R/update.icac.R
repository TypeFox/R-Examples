update.icac<-function(object,
		      what,
		      verbosity=2,
                      ...){
	S0<-as.data.frame(object$S0,stringsAsFactors=FALSE)
	S0$Trial<-as.numeric(as.character(object$data$Trial))
	if(verbosity>0)cat("updating ...\n")
	if(is.vector(what) && !is.list(what) && length(what)==3){
        		ic<-as.numeric(what[which(names(what)=="ic")])
        		trial<-as.numeric(what[which(names(what)=="trial")])
        		operation<-what[which(names(what)=="operation")]
		if(verbosity>1)cat(paste("  Trial = ", trial,
			"; IC = ",ic,"; operation = ",operation,
			" ...\n", sep=""))
            	if(operation=="+"){
            		newS <- as.data.frame(object$S,stringsAsFactors=FALSE)
            		newS$Trial <- as.numeric(as.character(object$data$Trial))
                	newS[newS$Trial==trial,ic]<-S0[S0$Trial==trial,ic]
            	}else{
            		newS <- as.data.frame(object$S,stringsAsFactors=FALSE)
                	newS$Trial <- as.numeric(as.character(object$data$Trial))
                	newS[newS$Trial==trial,ic]<-0
            	}
		if(verbosity>0)cat("re-mixing ICs ...")
		# re-mix corrected data = corrected source matrix * mixing matrix
            	S<-as.matrix(newS[,1:object$n.comp])
		tmp=as.data.frame(S%*%object$A)
		# give appropriate column names to corrected frame
		colnames(tmp)=object$channel 
		# replace EEG columns in original data frame with corrected ones
        		# and add channel mean (i.e., de-mean-center channels)
		for(ch in object$channel){ 
			if(verbosity>0)cat(".")
			object$data[,ch]=tmp[,ch]+object$col.means[ch]
		}
		if(verbosity>0)cat("\n")
        	}else if(is.list(what)){
            	S<-object$S
		for(kk in 1:length(what)){
			tmp.what<-what[[kk]]
        			ic<-as.numeric(tmp.what[1])
        			operation<-tmp.what[2]
			if(verbosity>1)cat(paste("  IC = ",ic,"; operation = ",
				operation," ...\n", sep=""))
            		if(operation=="+"){
                		S[,ic]<-S0[,ic]
            		}else{
                		S[,ic]<-0
            		}
		}
		if(verbosity>0)cat("re-mixing ICs ...")
		# re-mix corrected data = corrected source matrix * mixing matrix
		tmp<-as.data.frame(S%*%object$A)
		# give appropriate column names to corrected frame
		colnames(tmp)<-object$channel 
		# replace EEG columns in original data frame with corrected ones
        		# and add channel mean (i.e., de-mean-center channels)
		for(ch in object$channel){ 
			if(verbosity>0)cat(".")
			object$data[,ch]=tmp[,ch]+object$col.means[ch]
		}
		if(verbosity>0)cat("\n")
	}else{
		stop("argument \"what\" should be of a 3-element 
			 vector or a list of list of 2-element vectors\n")
	}
	if("updated"%in%names(object)){
		updated<-object$updated
		updated[[length(object$updated)+1]]<-what
	}else{
		updated<-list()
		updated[[1]]<-what
	}
	# create icac object
	res<-list(data=object$data,channel=object$channel,
		noise.sig=object$noise.sig,threshold=object$threshold,
		n.comp=object$n.comp,X=object$X,K=object$K,W=object$W,A=object$A,S=S,
		S0=object$S0,col.means=object$col.means,correlations=object$cor.info,
		correction.info=object$correction.info,updated=updated)
    if(!is.null(object$proctime)){
    	res$proctime<-object$proctime
    }

	class(res)<-"icac"

	# return output
	return(invisible(res))
}


