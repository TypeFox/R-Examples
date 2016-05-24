aug.x <- function(X,cp.locs,cp,delta = NULL){
X=as.matrix(X)
if(is.null(delta)) delta=cp.locs
if(length(cp.locs)!=length(delta) | length(cp.locs)!=sum(cp)){
	stop(paste("You must specify a correct changepoint structure!", 
		"\n"))
}
if(sum(cp)==0){
	X
} else{
new.x=NULL
sum.cp=cumsum(cp)
for(i in 1:length(cp)){
	X.i=matrix(X[,i],ncol=max(1,cp[i]),nrow=length(X[,1]))
	if(cp[i]!=0){
		CPs=cp.locs[(max(1,sum.cp[i-1]+1)):sum.cp[i]]
		SPs=delta[(max(1,sum.cp[i-1]+1)):sum.cp[i]]
		new.x=cbind(new.x,t(t(X.i)-SPs)*t(t(X.i)>CPs))
	}
}
cbind(X,new.x)
}
}

