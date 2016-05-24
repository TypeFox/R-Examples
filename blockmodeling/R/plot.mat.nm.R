"plot.mat.nm" <-
function(x=M,M=x,...,main.title=NULL,title.row="Row normalized",title.col="Column normalized",main.title.line=-2,par.set=list(mfrow=c(1,2))){
	if(is.null(main.title)){
		objName<-deparse(substitute(M))
		if(objName=="x")objName<-deparse(substitute(x))
		main.title <- paste("Matrix",objName)
	}
	if(!is.null(par)){
		par.def<-par(no.readonly = TRUE) 
		par(par.set)
	}
	row.normalized<-sweep(M, 1, apply(M, 1, sum),FUN="/")
	row.normalized[is.nan(row.normalized)]<-0
	plot.mat(M=row.normalized,main=title.row,outer.title=FALSE,...)
	column.normalized<-sweep(M, 2, apply(M, 2, sum),FUN="/")
	column.normalized[is.nan(column.normalized)]<-0
	plot.mat(M=column.normalized,main=title.col,outer.title=FALSE,...)
	title(main=main.title,outer=TRUE,line=main.title.line)
	if(!is.null(par.set))par(par.def)
}