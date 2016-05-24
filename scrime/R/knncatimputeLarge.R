`knncatimputeLarge` <-
function(data,mat.na=NULL,fac=NULL,fac.na=NULL,nn=3,distance=c("smc","cohen","snp1norm","pcc"),
		n.num=100,use.weights=TRUE,verbose=FALSE){
	if(is.null(mat.na)){
		rs<-rowSums(is.na(data))
		ids.na<-which(rs>0)
		if(length(ids.na)==0)
			stop("There are no missing values in data.")
		rn<-rownames(data)
		mat.na<-data[ids.na,,drop=FALSE]
		data<-data[-ids.na,,drop=FALSE]
		if(!is.null(fac)){
			fac.na<-fac[ids.na]
			fac<-fac[-ids.na]
		}
	}
	else{
		rs<-rowSums(is.na(mat.na))
		if(any(rs==0))
			stop("At least one of the rows of mat.na does not contain missing values.")
		rs<-rowSums(is.na(data))
		if(any(rs>0))
			stop("At least one of the rows of data contains missing values.")
		ids.na<-NULL
	}
	if(nn<1)
		stop("nn must be at least 1.")
	if(nn>nrow(data))
		stop("nn must be smaller than or equal to the number of rows of data.")
	if((is.null(fac) & !is.null(fac.na)) | (!is.null(fac) & is.null(fac.na)))
		stop("Either both or none of fac and fac.na has to be specified.")
	n.cat<-checkX1X2(data,mat.na)
	check4Monomorphism(data)
	check4Monomorphism(mat.na)
	if(is.null(fac)){
		fac<-rep(1,nrow(data))
		fac.na<-rep(1,nrow(mat.na))
		verbose<-FALSE
	}
	if(length(fac)!=nrow(data))
		stop("The length of fac must be equal to the number of columns of data.")
	if(length(fac.na)!=nrow(mat.na))
		stop("The length of fac.na must be equal to the number of columns of mat.na.")
	vec.split<-sort(unique(fac.na))
	tmp.split<-unique(fac)
	if(any(!vec.split%in%tmp.split))
		stop("At least one of the values in fac.na is not in fac.")
	distance<-match.arg(distance)
	distance<-paste(distance,"2Mats",sep="")
	for(i in 1:length(vec.split)){
		if(verbose)
			cat("Now considering factor ",vec.split[i],".",sep="")
		tmp.mat<-data[fac==vec.split[i],,drop=FALSE]
		tmp.matna<-mat.na[fac.na==vec.split[i],,drop=FALSE]
		out<-replaceNAs(tmp.mat,tmp.matna,nn=nn,distance=distance,n.num=n.num,use.weights=use.weights,
			n.cat=n.cat)
		mat.na[fac.na==vec.split[i],]<-out
		if(verbose)
			cat(" Done.\n",sep="")
	}
	if(is.null(ids.na))
		return(mat.na)
	mat<-matrix(0,length(rs),ncol(mat.na))
	mat[ids.na,]<-mat.na
	mat[-ids.na,]<-data
	rownames(mat)<-rn
	mat	
}

