binarize.array <- function(x,min.filter=NA,var.filter=0,fc.filter=0,na.filter=FALSE,log.base=NA,use.gap=FALSE){

	filter <- c()
	if(!is.na(min.filter)){
		cat(paste("filtering all rows with no values greater than",min.filter,"\n"))
		filter <- c(filter,which(apply(x,MARGIN=1,max,na.rm=TRUE)>min.filter))
	}
	if (var.filter>0){
		cat(paste("filtering ",var.filter*100,"% of rows with lowest variation \n",sep=""))
		sds <- apply(x,MARGIN=1,sd,na.rm=TRUE)
		sd.order <- sort(sds,decreasing=FALSE,index.return=TRUE)$ix
		filter <- c(filter,sd.order[1:floor(var.filter*nrow(x))])
        }
	if(fc.filter>0){
		cat(paste("filtering all rows with no fold-change greater than",fc.filter,"\n"))
		if(is.na(log.base)){
			fcs <- apply(x,MARGIN=1,function(y)max(y,na.rm=TRUE)/min(y,na.rm=TRUE))
			filter <- c(filter,which(fcs<fc.filter))
		}
		if(!is.na(log.base)){
			fcs <- apply(x,MARGIN=1,function(y)max(y,na.rm=TRUE)-min(y,na.rm=TRUE))
			filter <- c(filter,which(fcs<log(fc.filter,base=log.base)))
		}
	}
	if(na.filter){
		cat("filtering out all rows with missing values \n")
		filter <- c(filter,which(apply(x,MARGIN=1,function(y)sum(is.na(y)))>0))
	}

	unfiltered <- setdiff(1:nrow(x),filter)
	
	output <- array(0,dim=dim(x))
	cat(paste("applying cluster-based binarization to",length(unfiltered),"rows of data. This may take some time... \n"))
	if(use.gap) cat("using gap-statistic to determine cluster number. if this takes too long, try setting 'use.gap=FALSE' \n")
	output[unfiltered,] <- t(apply(x[unfiltered,],MARGIN=1,clusterDisc,use.gap=use.gap))
	rownames(output) <- rownames(x)
	colnames(output) <- colnames(x)
	output
}
