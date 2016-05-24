print.CPMCGLM <-function(x, ...)
{
	if(class(x)!="CPMCGLM") stop("The argument Object must be a class CPMCGLM ")
	cl <- x$call
	if (!is.null(cl)){
		cat("Call:\n")
		dput(cl)
		cat("\n")
	}	
	

	cat("Generalize Linear Model Summary \n")
	cat("	Family:","",x$family,"\n")
	cat("	Link:","",x$link,"\n")
	cat("	Number of subject:","",x$n,"\n")
	cat("	Number of adjustment variable:","",x$adj,"\n")
	cat("\n")

	cat("Resampling \n")
	cat("	N: ",x$N,"\n")

	n.pvalue<-round(x$naive.pvalue,4)
	if(is.numeric(x$exact.pvalue)==T){e.pvalue<-round(x$exact.pvalue,4)}
	else{e.pvalue<-x$exact.pvalue}
	b.pvalue<-round(x$bonferroni.adjusted.pvalue,4)
	pa.pvalue<-round(x$parametric.bootstrap.adjusted.pvalue,4)
	if(is.numeric(x$permutation.adjusted.pvalue)){pe.pvalue<-round(x$permutation.adjusted.pvalue,4)}
	else{pe.pvalue<-x$permutation.adjusted.pvalue}
	cat("\nBest coding \n")
	if(substr(x$BC,1,8)=="Quantile"){
		
		if (length(x$bestcod)==1){
		cat("Method:","","Dichotomous","","transformation","\n")}
		else{cat("Method:","",substr(x$BC,1,8),"\n")}
		cat("Value of the order quantile cutpoints:","",x$bestcod,"\n")
		cat("Value of the quantile cutpoints:","",x$vq,"\n")
	}
	if(substr(x$BC,1,8)=="Cutpoint"){
		cat("Method:","",substr(x$BC,1,8),"\n")
		cat("Value of the quantile cutpoints:","",x$bestcod,"\n")
	}
	if(substr(x$BC,1,6)=="Boxcox"){
		cat("Method:","",substr(x$BC,1,6),"\n")
		cat("Value of the BoxCox parameter:","",x$bestcod,"\n")
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

