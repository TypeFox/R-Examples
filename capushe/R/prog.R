AICcapushe = function(data,n){

if(!(is.numeric(n))){
	stop("n must be numeric")
	}

if (is.character(data)){data=read.table(data)}

if (length(data[1,])!=4){
	stop("data must have 4 columns with name, penshape, complexity and contrast")
	}
data=data.frame(data)
names(data)=c("model","pen","complexity","contrast")
if(any(is.na(data))){
	bad=which(is.na(data),arr.ind=TRUE)[,1]
	repbad=which(duplicated(bad),arr.ind=TRUE)
	if (length(repbad)!=0){bad=bad[-repbad]}
	data=data[-bad,]
	if (length(bad)==1){
		warning(paste("1 line has been removed by AICcapushe"))
		}
		else{warning(paste(as.character(length(bad)),"lines have been removed by AICcapushe"))}
	}
leng=length(data$model)
if (leng==0){
	stop("We need more correct observation")
	}
if ((length(data$pen)!=leng)|(length(data$complexity)!=leng)|(length(data$contrast)!=leng)){
	stop("The lengths of the columns must be equal")
	}
if (!prod(data$complexity>=0)){
	stop("Complexity must be positive")
	}
if (!(floor(n)==n)|(n<0)){
	stop("The number of observation must be an integer positive")
	}
if(n==Inf){
	stop("n must be different of Inf")
	}

AIC = data$contrast+data$complexity/n
modelmin=which.min(AIC)
result=list(data$model[modelmin],AIC[modelmin])
names(result)=c("model","AIC")
return(result)
}

BICcapushe = function(data,n){

if(!(is.numeric(n))){
	stop("n must be numeric")
	}

if (is.character(data)){data=read.table(data)}

if (length(data[1,])==2){
	data=matrix(c(data[,1],data[,1],data[,1],data[,2]),ncol=4)
	warning("The names of model are 'complexity'")
	}
if (length(data[1,])==3){
	data=matrix(c(data[,1],data[,2],data[,2],data[,3]),ncol=4)
	}
if (length(data[1,])!=4){
	stop("data must have 2, 3 or 4 columns")
	}
data=data.frame(data)
names(data)=c("model","pen","complexity","contrast")
if(any(is.na(data))){
	bad=which(is.na(data),arr.ind=TRUE)[,1]
	repbad=which(duplicated(bad),arr.ind=TRUE)
	if (length(repbad)!=0){bad=bad[-repbad]}
	data=data[-bad,]
	if (length(bad)==1){
		warning(paste("1 line has been removed by BICcapushe"))
		}
		else{warning(paste(as.character(length(bad)),"lines have been removed by BICcapushe"))}
	}
leng=length(data$model)
if (leng==0){
	stop("We need more correct observation")
	}
if ((length(data$pen)!=leng)|(length(data$complexity)!=leng)|(length(data$contrast)!=leng)){
	stop("The lengths of the columns must be equal")
	}
if (!prod(data$complexity>=0)){
	stop("Complexity must be positive")
	}
if (!(floor(n)==n)|(n<0)){
	stop("The number of observation must be an integer positive")
	}
if(n==Inf){
	stop("n must be different of Inf")
	}

BIC = data$contrast+data$complexity*log(n)/(2*n)
modelmin=which.min(BIC)
result=list(data$model[modelmin],BIC[modelmin])
names(result)=c("model","BIC")
return(result)
}



setGeneric("validation",function(x,data2,...){standardGeneric("validation")})