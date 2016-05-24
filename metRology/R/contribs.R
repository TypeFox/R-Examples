contribs<-function(object, scope, as.sd=FALSE, keep.sign=TRUE, simplify=TRUE, expand.dot=TRUE) {
	b.names<-row.names(object$budget)
	if( missing(scope) ) {
		snames <- scope <- b.names
	} else {
		if(class(scope)=="formula") {
			if("." %in% all.vars(scope) && expand.dot) {
				c.scope<-as.character(scope)
				c.scope[2] <- gsub("\\^ ?2", "", c.scope[2])
				c.scope[2] <- gsub("([+-]?)\\.([+-]?)", "\\1.^2\\2", c.scope[2])
				scope<-as.formula(c.scope)
			}
			snames<-fnames<-attr(terms(scope, data=object$cov),"term.labels")
			fnames <- unique(fnames, c(sub("(.*):(.*)","\\2:\\1",fnames)))
		} else if(class(scope)=="expression") {
			snames<-all.vars(scope)
		} else if(class(scope)=="character") {
			snames <- scope
		}
		snames <- snames[ snames %in% b.names ]
	}
	
	scope.index<-match(snames, b.names)
	x.names<-b.names[scope.index]
	ci<-object$budget$c[scope.index]
	covmat<-object$cov[scope.index, scope.index]
	cormat<-object$cor[scope.index, scope.index]
	contrib<-outer(ci,ci,"*") * covmat
	
	if(simplify) {
		indices<-which(abs(cormat)>2*.Machine$double.eps & 
					upper.tri(cormat), arr.ind=TRUE)
		x<-c(diag(contrib), 2*contrib[indices])
		names(x)<-c(x.names, paste(x.names[indices[,1]], 
				x.names[indices[,2]], sep=":"))
		if(class(scope)=="formula") x<-x[names(x) %in% fnames]
	} else {
		x<-contrib
	}

	if( as.sd ) {
		x.sd <- sqrt(abs(x))
		if(keep.sign) x.sd<-x.sd*sign(x) 
		return(x.sd)
	} else {
		return(x)
	}
}
