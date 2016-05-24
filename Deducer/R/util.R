perm<-function(vec,duplicates=FALSE) 
{
	if(!is.vector(vec))
		stop("vec must be a vector")
	n<-length(vec)
    if (n == 1) 
        return(as.matrix(vec[1]))
    else if (n < 2) 
        return(NULL)
    z <- matrix(1)
	if(all(!duplicated(vec)))
		duplicates=TRUE
    for (i in 2:n) {
        x <- cbind(z, i)
        a <- c(1:i, 1:(i - 1))
        z <- matrix(0, ncol = ncol(x), nrow = i * nrow(x))
        z[1:nrow(x), ] <- x
        for (j in 2:i - 1) {
            z[j * nrow(x) + 1:nrow(x), ] <- x[, a[1:i + j]]
        }
    }
    dimnames(z) <- NULL
	z<-apply(z,c(1,2),function(x) vec[x])
    if(!duplicates) z[!duplicated(z),] else z
}


table.to.data<-function(x){
	x<-as.matrix(x)
	d<-dim(x)
	res<-matrix(,ncol=2,nrow=sum(x))
	rn<-rownames(x)
	if(is.null(rn))
		rn<-1:d[1]
	cn<-colnames(x)
	if(is.null(cn))
		cn<-1:d[2]
	cnt<-1
	for(i in 1:d[1])
		for(j in 1:d[2]){
			if(x[i,j]==0)
				next
			tmp<-matrix(rep(c(rn[i],cn[j]),x[i,j]),ncol=2,byrow=TRUE)
			res[cnt:(cnt+x[i,j]-1),]<-tmp
			cnt<-cnt+x[i,j]
		}
	
	res<-as.data.frame(res)
	res[[1]]<-factor(res[[1]],levels=rn)
	res[[2]]<-factor(res[[2]],levels=cn)
	res
}
	


dich<-function(variables,data=NULL,cut=NULL,group1=NULL,group2=NULL){
	arguments <- as.list(match.call()[-1])
	variables<-eval(substitute(variables),data,parent.frame())
	if(length(dim(variables))<1.5){
		variables<-d(variables)
		fn<-arguments$variables
		names(variables)<-if(is.call(fn)) format(fn) else as.character(fn)
	}
	for(var in names(variables)){
		x<-variables[[var]]
		if(!is.null(cut)){
			cutpoint<-deparse(substitute(cut))
			x<- recode.variables(x , paste("Lo:",cutpoint," -> '<",cutpoint,"';",cutpoint,":Hi -> '>=",cutpoint,"';",sep=""))[[1]]
			x<-factor(x)
		}else{
			x<-factor(x)
			if(!is.null(group1) && is.null(group2)){
				name<-deparse(substitute(group1))
				x<-recode.variables(x, paste(name," -> '",name,"';else -> 'Not",name,"';"))[[1]]
			}
			if(is.null(group1) && !is.null(group2)){
				name<-deparse(substitute(group2))
				x<-recode.variables(x, paste(name," -> '",name,"';else -> 'Not",name,"';"))[[1]]
			}
			if(!is.null(group1) && !is.null(group2)){
				name1<-deparse(substitute(group1))
				name2<-deparse(substitute(group2))
				x<-recode.variables(x, paste(name1," -> '",name1,"';",name2," -> '",name2,
								"';else -> NA;"))[[1]]
			}
		}
		if(var == names(variables)[1]){
			result<-d(x)
			names(result)<-paste(var,".grp",sep="")
		}else
			result[[paste(var,".grp",sep="")]] <- x
	}
	result
}

d<-function(..., row.names = NULL, check.rows = FALSE,
                check.names = FALSE,
                stringsAsFactors = FALSE){
	data.frame(...,row.names=row.names,check.rows=check.rows,check.names=check.names,stringsAsFactors=stringsAsFactors)
}


get.objects<-function(cn=NULL,env = globalenv(),includeInherited=TRUE){
	objs <- ls(envir=env)
	if(is.null(cn))
		return(objs)
		
	l <- list()	
	for(obj in objs){
		if(!includeInherited)
			call <- paste("'",cn,"'"," == class(",obj,")[1]",sep="")
		else
			call <-paste("'",cn,"'"," %in% class(",obj,")",sep="")
		include <- try(eval(parse(text=call), envir=env),silent=TRUE)
		if(!inherits(include,"try-error") && include)
			l[[length(l)+1]] <- obj
	}
	result<- unlist(l)
	if(is.null(result))
		result <- character()
	result
}
