cor.matrix<-function(variables,with.variables,data=NULL,test=cor.test,...){
	arguments <- as.list(match.call()[-1])
	variables<-eval(substitute(variables),data,parent.frame())
	if(length(dim(variables))<1.5){
		variables<-d(variables)
		fn <- arguments$variables
		names(variables)<-if(is.call(fn)) format(fn) else as.character(fn)
	}
	if(missing(with.variables))
		with.variables <-variables
	else{
		with.variables<-eval(substitute(with.variables),data,parent.frame())
		if(length(dim(with.variables))<1.5){
			with.variables<-d(with.variables)
			fn <- arguments$with.variables
			names(with.variables)<-if(is.call(fn)) format(fn) else as.character(fn)
		}		
	}
	cors<-list()
	for(var1 in colnames(variables)){
		cors[[var1]]<-list()
		for(var2  in colnames(with.variables)){
			tmp<-na.omit(data.frame(as.numeric(variables[[var1]]),as.numeric(with.variables[[var2]])))
			names(tmp)<-c(var1,var2)
			cors[[var1]][[var2]]<-test(tmp[[1]],tmp[[2]],...)
			attr(cors[[var1]][[var2]],"N")<-nrow(tmp)
		}
	}
	class(cors)<-"cor.matrix"
	cors
}

print.cor.matrix<-function(x,digits=4,N=TRUE,CI=TRUE,stat=TRUE,p.value=TRUE,...){
	if(is.null(x[[1]][[1]]$conf.int))
		CI=FALSE
	n1<-length(x)
	n2<-length(x[[1]])
	label.width<-7
	num.rows<-6	
	result<-as.table(matrix(NA,nrow=n2*num.rows,ncol=n1))
	r.names<-names(x[[1]])
	r.name.width<-max(nchar(r.names))
	if(r.name.width>20){
		r.names<-formatC(r.names,width=20)
		r.name.width<-20
	}else
		r.names<-formatC(r.names,width=r.name.width)
	c.names<-names(x)
	if(max(nchar(c.names))>20)
		c.names<-format(names(x),width=15)
	colnames(result)<-c.names
	for(i in 1:n2)
		rownames(result)[i*num.rows-(num.rows-1)]<-paste(formatC(r.names[i],width=r.name.width),
				formatC("cor",width=label.width),sep="")
	
	
	rownames(result)[rep(1:num.rows,n2)>1.5]<-formatC(rep(c("N","CI*","stat**","p-value",
							paste(rep("-",r.name.width+label.width),collapse="")),n2),
			width=r.name.width+label.width)
	for(j in 1:n1){
		for(i in (1:n2)){
			result[i*num.rows-(num.rows-1),j]<-format(x[[j]][[i]]$estimate,digits=digits,...)
			if(!is.null(attr(x[[j]][[i]],"N")))
				result[i*num.rows-(num.rows-2),j]<-format(attr(x[[j]][[i]],"N"),...)
			if(names(x[[1]])[i]==names(x)[j])
				next
			if(!is.null(x[[j]][[i]]$conf.int))
				result[i*num.rows-(num.rows-3),j]<-paste("(",format(x[[j]][[i]]$conf.int[1],digits=digits,...),",",
						format(x[[j]][[i]]$conf.int[2],digits=digits,...),")",sep="")
			if(!is.null(x[[j]][[i]]$statistic)){
				if(!is.null(x[[j]][[i]]$parameter)){
					if(length(x[[j]][[i]]$parameter)==1)
						result[i*num.rows-(num.rows-4),j]<-paste(format(x[[j]][[i]]$statistic,digits=digits,...),
								" (",format(x[[j]][[i]]$parameter,digits=digits,...),")",
								sep="")
					else
						result[i*num.rows-(num.rows-4),j]<-paste(format(x[[j]][[i]]$statistic,digits=digits,...),
								" (",format(x[[j]][[i]]$parameter[1],digits=digits,...),",",
								format(x[[j]][[i]]$parameter[2],digits=digits,...),")",
								sep="")
				}else
					result[i*num.rows-(num.rows-4),j]<-format(x[[j]][[i]]$statistic,digits=digits,...)
			}
			if(!is.null(x[[j]][[i]]$p.value))
				result[i*num.rows-(num.rows-5),j]<-format(round(x[[j]][[i]]$p.value,digits),
						digits=digits,nsmall=digits,...)
		}
	}
	display<-if(CI) rep(TRUE,num.rows*n2) else rep(1:num.rows,n2)!=3
	display<-display & (if(N) rep(TRUE,num.rows*n2) else rep(1:num.rows,n2)!=2)
	display<-display & (if(stat) rep(TRUE,num.rows*n2) else rep(1:num.rows,n2)!=4)
	display<-display & (if(p.value) rep(TRUE,num.rows*n2) else rep(1:num.rows,n2)!=5)
	cat("\n",format(x[[1]][[1]]$method, width = getOption("width"), justify = "centre"), "\n\n")
	print(as.table(result[display,,drop=FALSE]))
	if(stat & !is.null(x[[1]][[1]]$statistic)){
		cat("\t**",names(x[[1]][[1]]$statistic)[1])
		if(!is.null(x[[1]][[1]]$parameter))
			if(length(x[[1]][[1]]$parameter)==1)
				cat(" (",names(x[[1]][[1]]$parameter[1]),")\n",sep="")
			else
				cat(" (",names(x[[1]][[1]]$parameter[1]),",",names(x[[1]][[1]]$parameter[2]),")\n",sep="")
		else
			cat("\n")
	}
	if(CI & !is.null(x[[1]][[1]]$conf.int))
		if(!is.null(attr(x[[1]][[1]]$conf.int,"conf.level")))
			cat("\t * ",attr(x[[1]][[1]]$conf.int,"conf.level")*100,"% percent interval\n\n",sep="")
		else
			cat("\n")
	if(!is.null(x[[1]][[1]]$alternative))
		cat("\tHA:",x[[1]][[1]]$alternative,"\n\n")
	else
		cat("\n")
}

as.matrix.cor.matrix<-function(x,...){
	n1<-length(x);
	n2<-length(x[[1]])
	mat<-matrix(NA,n2,n1)
	for(i in 1:n2)
		for(j in 1:n1)
			mat[i,j]<-x[[j]][[i]]$estimate
	colnames(mat)<-names(x)
	rownames(mat)<-names(x[[1]])
	mat
}