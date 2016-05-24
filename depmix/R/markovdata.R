
##################################
#                                #
#   MARKOVDATA OBJECT FUNCTION   #
#                                #
##################################

markovdata <- function(dat,itemtypes,nitems=length(itemtypes),ntimes=length(as.matrix(dat))/nitems,replicates=rep(1,length(ntimes)),inames=NULL,dname=NULL,xm=NA) {
	TYPES=c("continuous","categorical","count","covariate")
	itt = sapply(itemtypes,FUN=function(x){pmatch(tolower(as.character(x)),TYPES)})
	itemtypes[which(!is.na(itt))]=TYPES[itt[which(!is.na(itt))]]
	if(any(is.na(itemtypes))) stop("Itemtypes incorrectly specified.\n")
	if(is.data.frame(dat)) dat <- data.matrix(dat)
	if(!is.matrix(dat)) dat <- matrix(dat)
	if(!(dim(dat)[1]==sum(ntimes) && dim(dat)[2]==nitems)) dat <- matrix(unlist(dat),sum(ntimes),nitems,byrow=TRUE)
	if(!(length(dat)==(sum(ntimes)*nitems))) stop("Data length incompatible with ntimes and nitems.\n")
	dat=replace(dat,which(dat==xm),NA)
	# rearrange data such that items occur first, and then the covariates
	x=c(which(itemtypes!="covariate"),which(itemtypes=="covariate"))
	dat=dat[,x,drop=FALSE]
	itemtypes=itemtypes[x]
 	if(is.null(colnames(dat))) colnames(dat)=itemtypes
	if(!is.null(inames)) colnames(dat)=inames[x]
	if(is.null(dname)) dname=paste(nitems,"-item data",sep="")
	rownames(dat)=NULL
	attr(dat,"dname")=dname
	attr(dat,"itemtypes")=itemtypes
	attr(dat,"ntimes")=ntimes
	attr(dat,"replicates")=replicates
  	class(dat) = "md"
	return(dat)
}

# functions that return attributes of md objects
ntimes <- function(object) {return(attributes(object)$ntimes)}
itemtypes <- function(object) {return(attributes(object)$itemtypes)}
dname <- function(object) {return(attributes(object)$dname)}
replicates <- function(object){return(attributes(object)$replicates)}

# ... and functions of those same attributes
ncov <- function(object) {return(sum(as.logical(which(itemtypes(object)=="covariate"))))}
inames <- function(object) {return(colnames(object))}
nitems <- function(object) {return(dim(object)[2])}
ind <- function(object) {return(length(attributes(object)$ntimes))}

summary.md <- function(object, ...) {
	cat("Data set:                  ", dname(object), "\n")
	cat(" nr of items:              ", nitems(object), "\n")
	cat(" item type(s):             ", itemtypes(object), "\n")
	if(ncov(object)>0) 
	cat(" nr of covariates:         ", ncov(object), "\n")
	cat(" item name(s):             ", inames(object), "\n")
	if(ind(object)>5) 
	cat(" length(s) of series:      ", ntimes(object)[1:5], " ... \n")
	if(ind(object)<=5) 
	cat(" length(s) of series:      ", ntimes(object), " \n")
 	if(ind(object)>1) 
 	cat(" nr of independent series: ", ind(object), "\n")
#   	if(any(weights(object)!=1)) {
#  		if(ind(object)>5)
#  		cat(" case weights:             ", weights(object)[1:5]," ... \n")
#  		else 
#  		cat(" case weights:             ", weights(object)," \n")
#  	}	
   	cat(" data:              \n")
	print(object[1:3,])
}

print.md <- function(x, ...) {
	cat("Data set:                  ", dname(x), "\n")
	cat(" nr of items:              ", nitems(x), "\n")
	cat(" item type(s):             ", itemtypes(x), "\n")
	if(ncov(x)>0) 
	cat(" nr of covariates:         ", ncov(x), "\n")
	cat(" item name(s):             ", inames(x), "\n")
	if(ind(x)>5) 
	cat(" length(s) of series:      ", ntimes(x)[1:5], " ... \n")
	if(ind(x)<=5) 
	cat(" length(s) of series:      ", ntimes(x), " \n")
	if(ind(x)>1) 
	cat(" nr of independent series: ", ind(x), "\n")
# 	cat(" case weights:             ", head(weights(x)),"\n")
	cat(" data: \n")
	print(x[,], ...)
}

plot.md <- function(x, nitems=1:(min(5,dim(x)[2])),nind=1:(min(5,length(attributes(x)$ntimes))), ...) {
	dat=x[,nitems,drop=FALSE]
	ind<-length(attributes(x)$ntimes)
	if(ind==1) plot.ts(dat,plot.type="multiple",type="l",main=attributes(x)$dname, ...)
 	else {
		layout.show()
		# make dat into (multiple) timeseries objects
		ntimes=attributes(x)$ntimes
		dat2 <- list()
		for (i in nind) {
			if(i==1) bg=1 else bg=sum(ntimes[1:(i-1)])+1
			en=(sum(ntimes[1:i]))
			dat2[[i]] <- ts(dat[bg:en,,drop=FALSE])
		}
		y=dat2[[nind[1]]]
 		if(length(nind)>1) {
 			for(i in 2:length(nind)) y=cbind(y,dat2[[nind[i]]])
 		}
		colnames(y)=rep(colnames(dat),length(nind))
		plot.ts(y,plot.type="multiple",type="l",main=attributes(x)$dname, ...)
 	}
}
