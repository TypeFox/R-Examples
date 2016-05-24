two.sample.test<-function(formula,data=NULL,test=t.test,...){

	variables<-eval(formula[[2]],data,parent.frame())
	factor.var<-eval(formula[[3]],data,parent.frame())
	if(length(dim(variables))<1.5){
		variables<-d(variables)
		fn<-formula[[2]]
		names(variables)<-if(is.call(fn)) format(fn) else as.character(fn)
	}
	if(length(dim(factor.var))>1.5)
		factor.var <- factor.var[,1]
	
	x<-as.factor(factor.var)
	y<-variables
	if(length(levels(as.factor(as.integer(x))))!=2L) stop("factor must have 2 and only two levels") 
	tests<-list()
	tests[[1]]<-list(func=deparse(substitute(test)),level=levels(x))
	sample1<-!is.na(x) & x==levels(x)[1]
	sample2<-!is.na(x) & x==levels(x)[2]
	for(var.name in colnames(y)){
		tmp.x<-na.omit(y[sample1,var.name,drop=TRUE])
		tmp.y<-na.omit(y[sample2,var.name,drop=TRUE])
		if(length(tmp.y)<3L || length(tmp.x)<3L){
			warning(paste(var.name,"must have at least 3 observations per group. Dropping..."))
			next
		}
		tests[[var.name]]<-test(tmp.x,tmp.y,...)

	}
	result<-multi.test(tests)
	if(identical(t.test,test))
		colnames(result)[1:2]<-c(paste("mean of",tests[[1]]$level[1]),paste("mean of",tests[[1]]$level[2]))	
	result	
}

multi.test<-function(tests){
	combine.strings<-function(strings){
		if(length(strings)==1)
			return(as.character(strings))
		tmp.str<-"("
		for(i in 1:length(strings))
			tmp.str<-paste(tmp.str,strings[i],if(i<length(strings)) "," else "",sep="")
		tmp.str<-paste(tmp.str,")",sep="")
		tmp.str
	}
	n<-length(tests)
	result<-as.data.frame(matrix(NA,ncol=8,nrow=n-1))
	rownames(result)<-names(tests)[-1]
	colnames(result)<-c("estimate 1","estimate 2",
						"Difference",paste(attr(tests[[2]]$conf.int,"conf.level")*100,
											"% CI Lower",sep=""),
						paste(attr(tests[[2]]$conf.int,"conf.level")*100,
											"% CI Upper",sep=""),
						names(tests[[2]]$statistic),paste(combine.strings(names(tests[[2]]$parameter)),"",sep=""),"p-value")
	for(i in 2:n){
		if(!is.null(tests[[i]]$est)){
			if(length(tests[[i]]$est)==2){
				result[i-1,1]<-tests[[i]]$est[1]
				result[i-1,2]<-tests[[i]]$est[2]
				result[i-1,3]<-tests[[i]]$est[1]-tests[[i]]$est[2]
			}else{
				if(i==2) colnames(result)[1]<-names(tests[[i]]$est)
				result[i-1,1]<-tests[[i]]$est[1]
			}
		}
		if(!is.null(tests[[i]]$conf.int)){
			result[i-1,4]<-tests[[i]]$conf.int[1]
			result[i-1,5]<-tests[[i]]$conf.int[2]
		}
		result[i-1,6]<-tests[[i]]$statistic
		if(!is.null(tests[[i]]$parameter)){
			if(length(tests[[i]]$parameter)==1)
				result[i-1,7]<-tests[[i]]$parameter
			else{
				result[i-1,7]<-combine.strings(tests[[i]]$parameter)
			}
		}
		result[i-1,8]<-tests[[i]]$p.value
	}
	empty.cols<-apply(result,2,function(x)all(is.na(x)))
	result<-as.data.frame(result[,!empty.cols])
	attr(result,"method")<-tests[[2]]$method
	attr(result,"alternative")<-tests[[2]]$alternative
	if(!is.null(tests[[2]]$null.value)){
		attr(result,"null.value")<-tests[[2]]$null.value
	}
	class(result)<-c("multi.test","data.frame")
	result
}

print.multi.test<-function(x,...){
	cat(format(attr(x,"method"),width=getOption("width"),justify="centre"),"\n")
	class(x)<-"data.frame"
	print(x,...)
	if(is.character(attr(x,"alternative")))
		cat("  HA:",attr(x,"alternative"),"\n")
	if(!is.null(attr(x,"null.value")))
		cat("  H0: ",names(attr(x,"null.value"))[1],"=",attr(x,"null.value"),"\n")
}