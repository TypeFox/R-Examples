summary.CPMCGLM <-function(object, ...)
{
	if(class(object)!="CPMCGLM") stop("The argument Object must be a class CPMCGLM ")


	cat("Summary of CPMCGLM Package\n")
	n.pvalue<-round(object$naive.pvalue,4)
	if(is.numeric(object$exact.pvalue)==T){e.pvalue<-round(object$exact.pvalue,4)}
	else{e.pvalue<-object$exact.pvalue}
	b.pvalue<-round(object$bonferroni.adjusted.pvalue,4)
	pa.pvalue<-round(object$parametric.bootstrap.adjusted.pvalue,4)
	pe.pvalue<-round(object$permutation.adjusted.pvalue,4)

	cat("\nBest coding \n")
	if(substr(object$BC,1,8)=="Quantile"){
		cat("Method:","",substr(object$BC,1,8),"\n")
		cat("Value of the quantile cutpoints:","",object$bestcod,"\n")
	}
	if(substr(object$BC,1,8)=="Cutpoint"){
		cat("Method:","",substr(object$BC,1,8),"\n")
		cat("Value of the quantile cutpoints:","",object$bestcod,"\n")
	}
	if(substr(object$BC,1,6)=="Boxcox"){
		cat("Method:","",substr(object$BC,1,6),"\n")
		cat("Value of the BoxCox parameter:","",object$bestcod,"\n")
	}
	if(object$BC=="Original"){
		cat("The best association is found with the original continuous variable","\n")
	}
	
	cat("\n")
	cat("Corresponding adjusted pvalue: \n")
	cat("\n")
	if(is.numeric(e.pvalue)==T){
		R1<-data.frame(rbind(n.pvalue,e.pvalue,b.pvalue,pa.pvalue,pe.pvalue))
		R1[R1<0.0001]<-paste("<0.",paste(rep("0",3),collapse=""),"1",sep="")
		R1[R1>1]<- 1
		colnames(R1)<-rep("Adjusted pvalue",length(b.pvalue))
		rownames(R1)<-c("naive","exact","bonferroni","bootstrap","permutation")
		print(R1,digits=3)
	}else{
		R1<-data.frame(rbind(n.pvalue,b.pvalue,pa.pvalue,pe.pvalue))
		R1[R1==0]<-paste("<0.",paste(rep("0",3),collapse=""),"1",sep="")
		R1[R1>1]<- 1
		colnames(R1)<-rep("Adjusted pvalue",length(b.pvalue))
		rownames(R1)<-c("naive","bonferroni","bootstrap","permutation")
		print(R1,digits=3)
		cat("exact:","",e.pvalue)
		cat("\n")

	}	
}
