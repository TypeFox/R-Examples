treeperm<-function(x, ...) UseMethod("treeperm")

Permutation <- function(data,factor){
	factor<-as.factor(factor)
	p<-.Call("calculate_pvalue",data[factor==levels(factor)[1]],data[factor==levels (factor)[2]])
	result<-list(Statistics=p[3],Permutations=p[2],pvalue=p[1])
	class(result)<-"Permutation"
	result
}

KPermutation <- function(data,factor){
	factor<-as.factor(factor)
	p<-.Call("calculate_K_pvalue",data,as.integer(factor))
	result<-list(Statistics=p[3],Permutations=p[2],pvalue=p[1])
	class(result)<-"KPermutation"
	result
}

FStatistics<-function(data,factor){
	.Call("calculate_FStatistics",data,as.integer(factor))
}

ReducedStatistics<-function(data,factor){
	.Call("calculate_ReducedStatistics",data,as.integer(factor))
}

Ktreeperm<- function(data,factor,size){
	observed<-ReducedStatistics(data,factor)
	count<-0
	for(i in 1:size){
		s<-sample(data,length(data),replace=F)
		f<-ReducedStatistics(s,factor)
		if(f>=observed){
			count<-count+1
		}
	}
	result<-list(Statistics=FStatistics(data,factor),Permutations=size,pvalue=count/size,data=data,factor=factor)
	class(result)<-"KPermutation"
	result
}

#' @export
treeperm.default<-function(x,data,factor,type,size,...){
	if(!is.factor(factor)){
		warning('The type of factor must be factor')
		factor<-as.factor(factor)
	}
	if(!is.numeric(data)){
		stop('The type of data must be numeric')
	}
	if(length(levels(factor))==1){
		stop('You must have more than one group for permutation test')
	}
	if(type=="exact"){
		if(length(levels(factor))==2){
			result<-Permutation(data,factor)
		}else{
			result<-KPermutation(data,factor)
		}
	}else{
		if(type=="approximate"){
			result<-Ktreeperm(data,factor,size)
		}else{
			stop('You must specify parameter type')
		}
	}
	ran<-list(result=result,call=match.call(),data=data,factor=factor)
	class(ran)<-"treeperm"
	ran
}

#' @export
print.KPermutation<-function(x,...){
	cat("Total permutations\n")
	print(x$Permutations)
	cat("F statistics\n")
	print(x$Statistics)
	cat("P value\n")
	print(x$pvalue)
}

#' @export
print.Permutation<-function(x,...){
	cat("Total permutations\n")
	print(x$Permutations)
	cat("Observed statistics\n")
	print(x$Statistics)
	cat("P value\n")
	print(x$pvalue)
}

#' @export
print.treeperm<-function(x,...){
	print(x$result)
}

#' @export
summary.treeperm<-function(object,...){
	result<-list(call=object$call,information=object$result,data=object$data,factor=object$factor)
	class(result)<-"summary.treeperm"
	result
}

#' @export
print.summary.treeperm<-function(x,...){
	cat("Call:\n")
	print(x$call)
	cat("Group:\n")
	print(x$factor)
	cat("Data:\n")
	print(x$data)
	cat("Results:\n")
	print(x$information)
}

GetDistribution<-function(x,size){
	v<-rep(1:1,size)
	for(i in 1:size){
		s<-sample(x$data,length(x$data),replace=F)
		v[i]<-FStatistics(s,x$factor)
	}
	v
}

#' @export
plot.treeperm<-function(x,size,...){
	ran<-GetDistribution(x,size)
	r<-hist(ran,breaks=5,main="Permutations estimated by Monto carlo method",xlab="F statistics",ylab="Frequency",col = "lightblue")
	top<-max(r$counts)
	print(x$result$Statistics)
	points(FStatistics(x$data,x$factor),-top/150,type="p",pch=17,col="red",lwd=3)
	text(FStatistics(x$data,x$factor),-top/50,cex=.65,"Observed F statistics")
}
#' @export
treeperm.formula<-function(formula,frame=list(),type,size,...){
	mf<-model.frame(formula=formula,data=frame)
	data<-model.response(mf)
	factor<-as.factor(mf[,2])
	result<-treeperm.default(data=data,factor=factor,type=type,size=size,...)
	result$call<-match.call()
	result$formula<-formula
	result
}
