ts.bspline <-
function(expr,ts.point=NULL,data.predict=100,df=round(ncol(expr)*0.8)) {
	require(splines)
	gene.bs<-vector('list',length=nrow(expr))
	names(gene.bs)<-rownames(expr)
	
	if (is.null(ts.point)) ts.point=0:(ncol(expr)-1)
	ts.point.predict<-seq(min(ts.point),max(ts.point),length.out=data.predict)
	expr.predict<-matrix(0,nrow(expr),data.predict,dimnames=list(rownames(expr),as.character(ts.point.predict)))

	c<-0
	for (i in 1:length(gene.bs)) {
		if ((i*100)%/%length(gene.bs)>c) {
			cat(c,'% done.\n',sep='')
			c<-c+5
		}
		gene.bs[[i]]<-lm(expr~bs(timepoint,df=df),data=data.frame(timepoint=ts.point,expr=as.numeric(expr[i,])))
		expr.predict[i,]<-predict(gene.bs[[i]],data.frame(timepoint=ts.point.predict))
	}
	cat('\n')
	return(expr.predict)
}

