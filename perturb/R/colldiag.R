colldiag <- function(mod,scale=TRUE,center=FALSE,add.intercept=TRUE) {
	result <- NULL
	if (center) add.intercept<-FALSE
	if (is.matrix(mod)||is.data.frame(mod)) {
		X<-as.matrix(mod)
		nms<-colnames(mod)
	}
	else if (!is.null(mod$call$formula)) {
		X<-mod$model[,-1] # delete the dependent variable
	}
	X<-na.omit(X) # delete missing cases
	if (add.intercept) {
		X<-cbind(1,X) # add the intercept
		colnames(X)[1]<-"intercept"
	}
	X<-scale(X,scale=scale,center=center)

	svdX<-svd(X)
	svdX$d
	condindx<-svdX$d[1]/svdX$d

	Phi=svdX$v%*%diag(1/svdX$d)
	Phi<-t(Phi^2)
	pi<-prop.table(Phi,2)

	dim(condindx)<-c(length(condindx),1)
	colnames(condindx)<-"cond.index"
	rownames(condindx)<-1:nrow(condindx)
	colnames(pi)<-colnames(X)
	result$condindx<-condindx
	result$pi<-pi
	class(result)<-"colldiag"
	result
}

print.colldiag <- function(x,dec.places=3,fuzz=NULL,fuzzchar=".",...){
	stopifnot(fuzz>0 & fuzz<1)
	stopifnot(is.character(fuzzchar))
	stopifnot(nchar(fuzzchar)==1)
	fuzzchar<-paste(" ",fuzzchar,sep="")
	width<-dec.places+2
	pi<-formatC(x$pi,format="f",width=width,digits=dec.places)
	if (!is.null(fuzz )) {
		pi[pi < fuzz] <- fuzzchar
	}
	width<-max(nchar(trunc(max(x$condindx))))+dec.places+2
	condindx<-formatC(x$condindx,format="f",width=width,digits=dec.places)
	colnames(condindx)<-NULL
	cat("Condition\nIndex\tVariance Decomposition Proportions\n")
	print(noquote(cbind(condindx,pi)))
}
