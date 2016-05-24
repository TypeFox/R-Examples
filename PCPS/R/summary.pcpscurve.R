#' @rdname pcps.curve
#' @encoding UTF-8
#' @export
summary.pcpscurve<-function(object,probs=c(0.025,0.975),...){
    res<-list()
    res$call<-object$call
    res$curve.obs<-object$curve.obs
	if(!is.null(object$curve.null.ts)){
		N<-length(object$curve.null.ts)
		X<-matrix(NA,N,dim(object$curve.obs)[1])
		Y<-X
		for(i in 1:N){
			Ntemp<-length(object$curve.null.ts[[i]][,1])
			X[i,1:Ntemp]<-object$curve.null.ts[[i]][,1]
			Y[i,1:Ntemp]<-object$curve.null.ts[[i]][,2]
		}
		mean_X<-apply(X,2,mean,na.rm=T)
		mean_Y<-apply(Y,2,mean,na.rm=T)
		if(length(probs)!=2){
			stop("\n Only two values are accepted in probs \n")
		}		
		quantile_X<-apply(X,2,quantile,na.rm=T,probs=probs)
		quantile_Y<-apply(Y,2,quantile,na.rm=T,probs=probs)
		res_TS<-cbind(mean_X, mean_Y, quantile_X[1,],quantile_X[2,],quantile_Y[1,],quantile_Y[2,])
		rownames(res_TS)=rownames(object$curve.obs)
		colnames(res_TS)=c("Mean_Cum_PCPS_Eig","Mean_Coe_Det", paste(probs[1],"%_Conf_Cum_PCPS_Eig",sep=""), paste(probs[2],"%_Conf_Cum_PCPS_Eig",sep=""), paste(probs[1],"%_Conf_Coe_Det",sep=""), paste(probs[2],"%_Conf_Coe_Det",sep=""))
		res$null.model.ts<-res_TS
	}
	if(!is.null(object$curve.null.bm)){
		N<-length(object$curve.null.bm)
		X<-matrix(NA,N,dim(object$curve.obs)[1])
		Y<-X
		for(i in 1:N){
			Ntemp<-length(object$curve.null.bm[[i]][,1])
			X[i,1:Ntemp]<-object$curve.null.bm[[i]][,1]
			Y[i,1:Ntemp]<-object$curve.null.bm[[i]][,2]
		}
		mean_X<-apply(X,2,mean,na.rm=T)
		mean_Y<-apply(Y,2,mean,na.rm=T)
		if(length(probs)!=2){
			stop("\n Only two values are accepted in probs \n")
		}		
		quantile_X<-apply(X,2,quantile,na.rm=T,probs=probs)
		quantile_Y<-apply(Y,2,quantile,na.rm=T,probs=probs)
		res_BM<-cbind(mean_X, mean_Y, quantile_X[1,],quantile_X[2,],quantile_Y[1,],quantile_Y[2,])
		rownames(res_BM)=rownames(object$curve.obs)
		colnames(res_BM)=c("Mean_Cum_PCPS_Eig","Mean_Coe_Det", paste(probs[1],"%_Conf_Cum_PCPS_Eig",sep=""), paste(probs[2],"%_Conf_Cum_PCPS_Eig",sep=""), paste(probs[1],"%_Conf_Coe_Det",sep=""), paste(probs[2],"%_Conf_Coe_Det",sep=""))
		res$null.model.bm<-res_BM
	}
    return(res)
}