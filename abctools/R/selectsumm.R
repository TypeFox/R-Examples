selectsumm <-function (obs,param, sumstats, obspar=NULL, ssmethod=mincrit, verbose=TRUE,final.dens=FALSE, ...) {

if(!is.matrix(obs)|is.data.frame(obs)){
        obs<-matrix(obs,nrow=1)
}
if(!is.matrix(param)|is.data.frame(param)){
        param<-as.matrix(param)
}
if(!is.matrix(sumstats)|is.data.frame(sumstats)){
        sumstats<-as.matrix(sumstats)
}
if(!is.null(obspar)|is.data.frame(obspar)){ 
        if(!is.matrix(obspar)){ 
                obspar<-matrix(obspar,byrow=T,ncol=ncol(param))
        }
	if(nrow(obs)!=nrow(obspar)){
		stop("Please supply observed statistics and observed parameter matrices with the same number of rows!\n")
	}
}

if (!length(colnames(param))) {
        colnames(param) <- paste("P", 1:ncol(param), sep = "")
}
if (!length(colnames(sumstats))) {
        colnames(sumstats) <- paste("C", 1:ncol(sumstats), sep = "")
}

argl <- list(...)
sumsubs<-1:ncol(sumstats)
ssargind <-match(names(argl),"sumsubs")
ssargind <-which(!is.na(ssargind))
if (length(ssargind)>0){
	sumsubs<-eval(argl[[ssargind]])
}

sumstats<-sumstats[,sumsubs]
obs<-obs[,sumsubs]

if(!is.matrix(obs)){
        obs<-matrix(obs,nrow=1)
}
if(!is.matrix(sumstats)){
        sumstats<-as.matrix(sumstats)
}

ndatasets <- nrow(obs)
nstats<-ncol(obs)
   
err <- critvals <- best<-vals<-NULL

record<-matrix(0,ndatasets,length(sumsubs))
colnames(record)<-paste("S",sumsubs,sep="")
rownames(record)<-rep("",ndatasets)

	for (i in 1:ndatasets) {
	#	if(verbose){
        		cat("dataset...", i, "\n")
	#	}
		resi<-ssmethod(obs[i,],param, sumstats, obspar[i,], verbose = verbose,final.dens=final.dens, ...)
		if(!is.null(resi$best)){
			record[i,resi$best]<-1
			rownames(record)[i]<-rownames(resi$best)
		}
        	if(!is.null(resi$err)){    
    		    err <- rbind(err, resi$err)
        	}
        	if(!is.null(resi$critvals)){    
        		critvals <- rbind(critvals, resi$critvals)
        	}
        	if (final.dens) {
			if(ncol(resi$post.sample)==1){
                		vals <- cbind(vals, resi$post.sample)
			}
			else{
                		vals <- abind(vals, resi$post.sample,along=3)
			}
        	}
	}

l<-list()

if (length(ssargind)!=0){
	l$sumsubs<-sumsubs
}

if(!is.null(critvals)){    
        l$critvals<-critvals
}
if(!is.null(err)){    
        l$err<-err
}
if(final.dens){
        l$post.sample<-vals
}

if(!is.null(resi$order)){
	l$posssubs<-resi$posssubs
}
if(!is.null(resi$sainfo)){
	l$sainfo<-resi$sainfo
}

if(sum(record)!=0){
	l$best<-record
}

return(l)

}

