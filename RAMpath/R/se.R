ramEffectSE<-function(object, effect, path=TRUE){	
	l.object<-object$ram$lavaan
	vcov.est<-vcov(l.object)
	## parameter table
	parTable<-l.object@ParTable
	parEst<-l.object@Fit@est

	ovname<-l.object@Data@ov$name
	
	ramEffectDirect<-function(v1, v2, path=TRUE, ovname, parTable, parEst, vcov.est){
	## the standard error is readily available
		if (path){
		## For path parameter
			if (v1 %in% ovname){
				## the path is a regression
				rname<-paste(v2, "~", v1, sep="")
				pos<-which(parTable$lhs == v2 & parTable$rhs == v1)
				if (parTable$label[pos]!="") rname<-parTable$label[pos]
		
				## get the position in vcov
				pos.cov<-which(rownames(vcov.est)==rname)
				effect.var<-vcov.est[pos.cov, pos.cov]
				est<-parEst[parTable$id[pos]]
			}else{
				## the path is a factor loading
				rname<-paste(v1, "=~", v2, sep="")
				pos<-which(parTable$lhs == v1 & parTable$rhs == v2)
				if (parTable$label[pos]!="") rname<-parTable$label[pos]
		
				## get the position in vcov
				pos.cov<-which(rownames(vcov.est)==rname)
				effect.var<-vcov.est[pos.cov, pos.cov]
				est<-parEst[parTable$id[pos]]
			}
		}else{
			## for variance / covariance parameters		
			pos<-which((parTable$lhs == v2 & parTable$rhs == v1)|(parTable$lhs == v1 & parTable$rhs == v2))
			rname<-paste(parTable$lhs[pos], "~~", parTable$rhs[pos], sep="")
			if (parTable$label[pos]!="") rname<-parTable$label[pos]		
			## get the position in vcov
			pos.cov<-which(rownames(vcov.est)==rname)
			effect.var<-vcov.est[pos.cov, pos.cov]
			est<-parEst[parTable$id[pos]]
		}
		return(c(est, effect.var, pos.cov))
	}
	
	if (path){
		effect.sub<-unlist(strsplit(effect, " > ", fixed=TRUE))
		effect.len<-length(effect.sub)

		if (effect.len<=2){
			## the standard error is readily available
			v1<-effect.sub[1]
			v2<-effect.sub[2]	
			effect.var<-ramEffectDirect(effect.sub[1], effect.sub[2], path=path, ovname, parTable, parEst, vcov.est)
			effect.se<-sqrt(effect.var[2])
		}else{
			temp.name<-temp.est<-temp.pos<-NULL
			ind.name<-''
			for (i in 1:(effect.len-1)){
				temp.effect<-ramEffectDirect(effect.sub[i], effect.sub[i+1], path=path, ovname, parTable, parEst, vcov.est)
				temp.name<-c(temp.name, paste('a', i, sep=''))
				ind.name<-paste(ind.name, 'a', i, '*', sep='')
				temp.est<-c(temp.est, temp.effect[2])
				temp.pos<-c(temp.pos, temp.effect[3])		
			}
			ind.name<-substr(ind.name, 1, nchar(ind.name)-1)
	
			## now take the first derivative
			ind.exp<-parse(text=ind.name)
			names(temp.est)<-temp.name
			par.list<-as.list(temp.est)	

			ind.deriv<-deriv(ind.exp, temp.name)
			first.deriv<-eval(ind.deriv, par.list)
			first.deriv<-attributes(first.deriv)$gradient
	
			var.par<-vcov.est[temp.pos, temp.pos]
			var.ind<-first.deriv%*%var.par%*%t(first.deriv)
			effect.se<-sqrt(var.ind)
		}
	}
	effect.se
}