`predict.pamCat` <-
function(object,newdata,theta=NULL,add.nvar=FALSE,type=c("class","prob"),...){
	if(missing(newdata))
		stop("newdata must be specified.")
	if(any(is.na(newdata)))
		stop("No missing values allowed in newdata.")
	n.cat<-object$n.cat
	if(any(!newdata%in%1:n.cat))
		stop("newdata must consist of values between 1 and ",n.cat,".") 
	mat.chisq<-object$mat.chisq
	if(nrow(mat.chisq)!=nrow(newdata))
		stop("The number of rows of newdata differs from the number of variables.")
	if(is.null(rownames(newdata)))
		stop("newdata must have row names.")
	if(any(rownames(mat.chisq)!=rownames(newdata)))
		stop("Each of the rows of newdata must represent the same variable\n",
			"that is contained in the corresponding row of data.")
	if(is.null(theta)){
		theta<-rev(object$mat.theta[,1])[which.min(rev(object$mat.theta[,3]))]
		warning("Since theta has not been specified, it is set to ",theta,".")
	}
	if(theta<=0)
		stop("theta must be stricly positive.")
	if(theta>=max(mat.chisq)){
		pred<-rep(NA,ncol(newdata))
		out<-if(add.nvar) list(pred=pred,n.var=0) else pred
		warning("Predicted values are set to NA, since theta is larger than\n",
			"the maximum value of the test statistics.")
		return(out)
	}
	new.stats<-mat.chisq-theta
	new.stats[new.stats<0]<-0
	ids<-rowSums(new.stats)>0
	nu<-sqrt(new.stats[ids,]/mat.chisq[ids,])
	N<-object$mat.obs[ids,,drop=FALSE]
	Nexp<-object$mat.exp[ids,,drop=FALSE]
	Ntheta<-Nexp+as.vector(nu)*(N-Nexp)
	nr.obs<-as.numeric(object$tab.cl)
	prior<-nr.obs/sum(nr.obs)
	n.lev<-length(nr.obs)
	n.var<-sum(ids)
	if(length(unique(nr.obs))==1)
		nr.obs<-unique(nr.obs)
	else
		nr.obs<-rep(nr.obs,e=n.var)
	Ntheta<-Ntheta/nr.obs
	n.new<-ncol(newdata)
	type<-match.arg(type)
	pred<-if(type=="class") numeric(n.new) else matrix(0,n.new,n.lev)
	newdata<-newdata[ids,,drop=FALSE]
	mat.prob<-matrix(0,n.var,n.lev)
	for(i in 1:n.new){
		for(j in 1:n.cat){
			tmp.ids<-which(newdata[,i]==j)
			if(length(tmp.ids)>0)
				mat.prob[tmp.ids,]<-Ntheta[tmp.ids,n.lev*(j-1)+(1:n.lev)]
		}
		prob<-apply(mat.prob,2,prod)
		prob<-prob*prior
		if(type=="class")
			pred[i]<-which.max(prob)
		else
			pred[i,]<-prob/sum(prob)
	}
	if(type=="prob")
		colnames(pred)<-1:n.lev
	if(add.nvar)
		return(list(pred=pred,n.var=n.var))
	pred
}

