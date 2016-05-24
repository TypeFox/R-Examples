
# Corey M. Yanofsky, February 19, 2009
# M. Padilla, (changes) 2013 
# x is treatment
# y is control
#---------------
sameAsX_names<-function(x=NULL,y=NULL){
	
	assert.are(list(x,y),c("numeric","logical","character"))
	#stopifnot(length(x)>=length(y) && all(names(y)%in%names(x)))
	z<-sameXY_names(x=x,y=y)
	
	z$ny[names(x)]
}
sameAsY<-function(x=NULL,y=NULL){sameAsX_names(x=y,y=x)}
k.matrix2prnSet<-function(x,y){
	if(!is_any(x,c("matrix","nxprnSet"))||!is_any(y,c("matrix","NULL","nxprnSet"))){stop("bad input class")}
	if(is(x,"matrix")){x<-new_nxprnSet(exprs =x)}
	if(is(y,"matrix")){y<-new_nxprnSet(exprs =y)}
	assert.is(x,"nxprnSet");assert.is(y,c("nxprnSet","NULL"))
	list(x=x,y=y)
}
get_testfun<-function(x,y=NULL,test.fun=t.test_safe,opt=c("statistic","p.value"),...){opt<-opt[1]
	z<-k.xy.prnset2matrix(x=x,y=y,paired=F)
	x<-z$x;y<-z$y
	if(is(test.fun,'character')){test.fun<-eval(parse(text=test.fun))}

	fct<-function(i){z<-as.numeric(NA)
		if(length(y)>0){z<-try(test.fun(x=x[i,],y=y[i,],...))}
		else{z<-try(test.fun(x=x[i,],...))}
		if(is_err(z)){z<-as.numeric(NA)}
		else{z<-z[[opt]]}
		z	
	}

	zo<-vapply(1:nrow(x),FUN.VALUE=numeric(1),FUN=fct)
	names(zo)<-rownames(x)
	zo
}
get_pvalues<-function(x,y=NULL,test.fun=t.test_safe,...){
	get_testfun(x=x,y=y,test.fun=test.fun,opt=c("p.value"),...)
}
get_stats<-function(x,y=NULL,test.fun=t.test_safe,...){
	get_testfun(x=x,y=y,test.fun=test.fun,opt=c("statistic"),...)
}
get_pvalfun<-function(pval.fun=t.test,alternative='greater',arglis.pvalfun=list(),...){
	
	zx<-paste("function(x,y=NULL,alternative='",alternative,"'){do.call(pval.fun,c(list(x=x,y=y,alternative=alternative),arglis.pvalfun))}",sep="")
	pfx <- eval(parse(text =zx))
	if(is(pfx,'character')){pfx<-eval(parse(text=paste(pfx,sep="")))}
	stopifnot(is(pfx,'function'))
	PValueFUN(FUN=pfx,alternative=alternative,...)
}

interest.param.fmeans<-function(x,y=NULL,paired=F,param=NULL,param0=NULL,logx=TRUE){
	z<-nprnSet2matrix(x=x,y=y,paired=F)
	x<-z$x;y<-z$y
		
	xbar <- rowMeans(x,na.rm = TRUE)
	if(is(y,'matrix')){ybar <- rowMeans(y,na.rm = TRUE)}
	else{ybar <- 0}

	if(logx&&length(param)==0){param<-xbar-ybar}
	if(logx&&length(param0)==0){param0<-0}
	if(!logx&&length(param)==0){param<-xbar/ybar}
	if(!logx&&length(param0)==0){param0<-1}
	param<-MakeNames(param)
	param<-sameAsY(x=param,y=xbar)

list(xbar=xbar,ybar=ybar,param=param,param0=param0)
	
}
k.checkings.locfdr<-function(pvalue=NULL,stat=NULL,param=NULL,nulltype){
	assert.is(nulltype,'numeric')
	stopifnot(nulltype%in%c(0:3))
	stopifnot(length(stat)>0||length(pvalue)>0)

	if(length(pvalue)>0){
		stopifnot(all(is_prob(pvalue[is.finite(pvalue)])))
		stat<-qnorm(pvalue)
		pvalue<-MakeNames(pvalue)}
	if(is.null(param)){warning("empty param");param<-rep(as.numeric(NA),length(stat))}
	stat<-MakeNames(stat)
	param<-MakeNames(param)
	
	orig.pvalue<-pvalue
	orig.stat<-stat

	if(length(pvalue)>0){stopifnot(length(pvalue)==length(stat) && identical(names(pvalue),names(stat)))}
	stopifnot(length(param)==length(stat) && identical(names(param),names(stat)))
	
	ind<-which(is.finite(stat))
	stat<-stat[ind]
	if(length(pvalue)>0){pvalue<-pvalue[ind]}
	#param<-param[ind]
	
	list(pvalue=pvalue,stat=stat,param=param,orig.pvalue=orig.pvalue,orig.stat=orig.stat,nulltype=nulltype)
}



#----------------
t.test_safe <- function(x,...){
	junk<- try(t.test(x=x,...), silent = TRUE)
	if(is(junk, "try-error") && substr(junk,nchar(junk) - 29, nchar(junk) - 1) == "data are essentially constant")
	{	
		junk <- list()
		junk$p.value <- 0
	}
	if(is.na(junk$p.value)) junk$p.value <- 1
	junk
}

w.test_safe <- function(x,...){
	junk<- try(wilcox.test(x=x,...), silent = TRUE)
	if(is(junk, "try-error") && substr(junk,nchar(junk) - 29, nchar(junk) - 1) == "data are essentially constant")
	{	
		junk <- list()
		junk$p.value <- 0
	}
	if(is.na(junk$p.value)) junk$p.value <- 1
	junk
}

locfdr.really_safe <- function(...){failed.num<-as.numeric(NA)#old:failed.num<-0
	ctr = 1;arglis<-formalArgs(def=locfdr)
	temp <- list(...);nms<-names(temp)
	ix<-c(which(nms=="zz"),1)[1]
	ozz<-temp[[ix]]
	
	temp$plot<-0
	temp<-temp[names(temp)%in%arglis]
	out <- try(do.call(locfdr,temp)$fdr, silent = TRUE)

	while(ctr < 10 && class(out) == "try-error")
	{
		temp <- list(...);nms<-names(temp)
		if(length(nms)==0||(length(nms)>0&&nchar(nms[1]) == 0))#if(names(temp)[1] == "")#old:if(names(temp)[[1]] == "")
		{
			junk <- names(temp)
			junk[1] = "zz"
			names(temp) <- junk
		}
		ctr = ctr + 1
		idx <- sort.list(abs(temp$zz), decreasing = TRUE)[1:10*ctr]
		idx2 <- setdiff(1:length(temp$zz), idx)
		out <- rep(failed.num,length(temp$zz))
		temp$zz <- temp$zz[idx2]
		fdr <- try(do.call(locfdr.really_safe,c(temp)),silent=T)
		#fdr <- try(eval(as.call(c(locfdr.really_safe,temp))),silent=T)
		
		out[idx2] <- fdr
		out[idx] <- failed.num
	}
	if(!is(out,'numeric')){stop('locfdr failed')}
	names(out)<-names(ozz);out
}
	

locfdr.safe <- function(...){failed.num<-as.numeric(NA)#old:failed.num<-0
	
	ctr = 1;arglis<-formalArgs(def=locfdr)
	temp <- list(...);nms<-names(temp)
	ix<-c(which(nms=="zz"),1)[1]
	ozz<-temp[[ix]]
	
	temp$plot<-0
	temp<-temp[names(temp)%in%arglis]
	out <- try(do.call(locfdr,temp)$fdr, silent = TRUE)
	
	if(length(nms)==0||(length(nms)>0&&nchar(nms[1]) == 0))
	{
		junk <- names(temp)
		junk[1] = "zz"
		names(temp) <- junk
	}
	idx <- which(!is.finite(temp$zz))#old:is.infinite(temp$zz)
	if(length(idx) > 0)
	{
		idx2 <- setdiff(1:length(temp$zz), idx)
		out <- rep(failed.num,length(temp$zz))
		temp$zz <- temp$zz[idx2]
		fdr <- eval(as.call(c(locfdr.really_safe,temp)))
		out[idx2] <- fdr
		out[idx] <- failed.num
	}
	else out <- locfdr.really_safe(...)
	names(out)<-names(ozz);out
}
 
estimator_func_generator <- function(...){
	function(x, y,...){
		if(missing(y)||is.null(y)){p <- Prob0(x = x, ...)}
		else{p <- Prob0(x = x, y = y, normalization = NULL, ...)}
	
		(1 - p@.Data)*p@Location1
	}
}
#---------------
locfdr_empnull_w_test_estimator <- estimator_func_generator(probabilityNull = eFdr, PValue = PValueFUN(wilcox.test))
locfdr_empnull_t_test_estimator <- estimator_func_generator(probabilityNull = eFdr, PValue = PValueFUN(t.test))
locfdr_theornull_w_test_estimator <- estimator_func_generator(probabilityNull = eFdr0, PValue = PValueFUN(wilcox.test))
locfdr_theornull_t_test_estimator <- estimator_func_generator(probabilityNull = eFdr0, PValue = PValueFUN(t.test))
locfdr_moderated_t_estimator <- estimator_func_generator(probabilityNull = eFdr0, PValue = sFdr)
limma_estimator <- estimator_func_generator(probabilityNull = sFdr)
pseudo_estimator <- estimator_func_generator(probabilityNull = pseudoFdr)
threshold0.5_estimator <- estimator_func_generator(probabilityNull = pseudoFdr, foldchange.threshold = 0.5)
threshold1_estimator <- estimator_func_generator(probabilityNull = pseudoFdr, foldchange.threshold = 1)
threshold2_estimator <- estimator_func_generator(probabilityNull = pseudoFdr,foldchange.threshold = 2)
#
scottbergermodel_MCMC <- function(obj, n_iter = 10, burnin = 2, nchains = 15, opt_init_val){
	xbar <- obj[["xbar"]]
	sample_size_factor <- obj[["sample_size_factor"]]
	logtarget <- function(x) scottbergermodel(xbar, sample_size_factor, x[1], x[2], x[3])
	if(missing(opt_init_val)) opt_init_val <- c(log(var(xbar)/20),log(var(xbar)/20),log(9))
	ctr <- 19
	while(is.infinite(logtarget(opt_init_val)) && ctr > 0)
	{
		ctr <- ctr - 1
		opt_init_val <- c(log(var(xbar)/ctr),log(var(xbar)/ctr),log(9))
	}
	optpoint <- optim(par=opt_init_val, logtarget, hessian = TRUE, control = list(fnscale = -1))
	init_dist_chol <- chol(-solve(optpoint$hessian))
	n_param <- 3
	init_state <- vector("list")
	for(i in 1:nchains) init_state[[i]] <- as.vector(optpoint$par + t(t(init_dist_chol)%*%matrix(rnorm(n_param),nrow = n_param)))
	state <- DEMC_sample(logtarget,burnin + n_iter,init_state)
	obj <- scottbergermodel_posterior_mean(xbar, sample_size_factor, state[(burnin +(1:n_iter))])
	(1 - obj[["p0"]])*obj[["mu_given_DE"]]
}

scottbergermodel_estimator <- function(x, y=NULL, ...){
	scottbergermodel_MCMC(obj = get_xbar_and_sample_size_factor(x=x,y=y), ...)
	}

scottbergermodel_scaled_data_estimator <- function(x, y=NULL, n_iter=10, burnin=2, nchains=15, ...){
	
	z<-k.xy.prnset2matrix(x=x,y=y,paired=F)
	x<-z$x;y<-z$y
		
	if(length(y)>0){
		
		if(nrow(x)!=nrow(y)) stop ("number of features (rows) are not same")
		
		scale_factor <- sapply(1:nrow(x), function(i) {
			dat1 <- na.omit(as.numeric(x[i,]))
			dat2 <- na.omit(as.numeric(y[i,]))
			n1 <- length(dat1)
			n2 <- length(dat2)
			((var(dat1)*(n1-1)) + (var(dat2)*(n2 - 1)))/(n1+n2-2)})
		
	}
	else{
		scale_factor <- sapply(1:nrow(x), function(i) var(na.omit(as.numeric(x[i,]))))
	}
	#----------------------
	obj <- get_xbar_and_sample_size_factor(x=x,y=y)
	obj$scale_factor <- sqrt(scale_factor + median(scale_factor,na.rm = TRUE))
	obj$xbar <- obj$xbar/obj$scale_factor
 	states <- scottbergermodel_posterior_sample(obj, n_iter, burnin, nchains, ...)
 	scottbergermodel_scaled_data_posterior_mean(xbar=obj$xbar, sample_size_factor=obj$sample_size_factor,
						    scale_factor=obj$scale_factor, states=states)

	
	}


##-----------------
old.locfdr_estimator <- function(x, y=NULL,pval, p0,...){
	k.pval<-function(x,y=NULL){
		if(length(y)==0){
			k.fun<-'try(t.test_safe(x = x[i,], alternative = "less")$p.value))'
		}
		else{
			k.fun<-'try(t.test_safe(x = x[i,], y = y[i,], alternative = "less")$p.value))'
		}
		pval <- vapply(1:nrow(x),FUN.VALUE=numeric(1),FUN= function(i){
			eval(parse(text=paste("zaux<-",k.fun,sep="")))
			if(is_err(zaux)){zaux<-as.numeric(NA)}
			zaux})
		pval
	}
	#-----
	z<-k.xy.prnset2matrix(x=x,y=y,paired=F)
	x<-z$x;y<-z$y
	xbar <- rowMeans(x,na.rm = TRUE)
	
	if(length(y)>0){
		ybar <- rowMeans(y,na.rm = TRUE)
	}
	else{ybar<-0}
	
	if(missing(pval) && missing(p0)){pval<-k.pval(x=x,y=y)}
	if(missing(p0)) {p0 <- locfdr.safe(qnorm(pval), nulltype = 0,...)}
		(1 - p0)*(xbar - ybar)
	
}

k.est.sev<-function(x, y=NULL,opt="basic_BF_estimator", ...){opt<-opt[1]
	
	k.x.fun<-function(x,fun){
		assert.is(fun,"character");assert.is(x,"matrix")
		kfct<-function(mat,...){fun(as.numeric(mat[1:ncol(x)]),...)}#{eval(parse(text=paste(fun,"as.numeric(mat[1:ncol(x)]),...))",sep="")))}
		p0 <- apply(x, 1, FUN=kfct)
		p0
	}
	k.xy.fun<-function(x,y,fun){
		assert.is(fun,"character");assert.is(x,"matrix");assert.is(y,"matrix")
		kfct<-function(mat,...){fun(as.numeric(mat[1:ncol(x)]),as.numeric(mat[((1:ncol(y))+ncol(x))]),...)}
		#eval(parse(text=paste(fun,"as.numeric(mat[1:ncol(x)]),as.numeric(mat[((1:ncol(y))+ncol(x))]),...))",sep="")))
		p0 <- apply(cbind(x,y), 1, FUN=kfct)
		p0
	}
	#---
	if(opt%in%"basic_BF_estimator"){k.fun<-"logBF2posprob(BF.one.covariate("}
	if(opt%in%"AICc_estimator"){k.fun<-"logBF2posprob(.5*AICc.one.covariate("}
	if(opt%in%"BIC_estimator"){k.fun<-"logBF2posprob(.5*BIC.one.covariate("}
		
	#-----
	z<-k.xy.prnset2matrix(x=x,y=y,paired=F)
	x<-z$x;y<-z$y
	xbar <- rowMeans(x,na.rm = TRUE)
	
	if(length(y)>0){
		ybar <- rowMeans(y,na.rm = TRUE)
		p0 <- k.xy.fun(x=x,y=y,fun=k.fun)
	}
	else{
		ybar<-0
		p0 <- k.x.fun(x=x,fun=k.fun)
	}		
	(1 - p0)*(xbar - ybar)
}

basic_BF_estimator <- function(x, y=NULL, ...){k.est.sev(x=x, y=y,opt="basic_BF_estimator", ...)}	
AICc_estimator <- function(x, y=NULL, ...){k.est.sev(x=x, y=y,opt="AICc_estimator", ...)}
BIC_estimator <- function(x, y=NULL, ...){k.est.sev(x=x, y=y,opt="BIC_estimator", ...)}
old.EApvalC<- function(x,y=NULL,method="LFDR",pval, alpha=0.05,...){
	k.pval<-function(x,y=NULL,alternative = "less"){
		if(length(y)==0){
			k.fun<-'try(t.test_safe(x = x[i,], alternative = alternative)$p.value))'
		}
		else{
			k.fun<-'try(t.test_safe(x = x[i,], y = y[i,], alternative = alternative)$p.value))'
		}
		pval <- vapply(1:nrow(x),FUN.VALUE=numeric(1),FUN= function(i){
			eval(parse(text=paste("zaux<-",k.fun,sep="")))
			if(is_err(zaux)){zaux<-as.numeric(NA)}
			zaux})
		pval
	}
	#-----
	z<-k.xy.prnset2matrix(x=x,y=y,paired=F)
	x<-z$x;y<-z$y
	xbar <- rowMeans(x,na.rm = TRUE)
	
	if(length(y)>0){
		ybar <- rowMeans(y,na.rm = TRUE)
	}
	else{ybar<-0}
	
	if(missing(pval)&&method %in% c("TWER","FWER")){pval<-k.pval(x=x,y=y,alternative="two.sided")}
	else if(missing(pval)&&method %in% c("LFDR")){pval<-k.pval(x=x,y=y,alternative="less")}
	else if(missing(pval)&&!method %in% c("LFDR","TWER","FWER")){stop("bad method")}

	stopifnot(length(pval[is.finite(pval)]>0))	

	
	if(method == "FWER"){pval <- mt.rawp2adjp(pval, proc = "SidakSD")$adjp[,2]}
	else if(method == "LFDR"){pval <- locfdr.safe(qnorm(pval), nulltype = 0)}
	else{stop("bad method")}

	
	warn_state <- getOption("warn")
	options(warn = -1)
	if(length(y)==0){xbar[pval > alpha] <- 0;ybar<-0}
	else{
		xbar[pval > alpha] <- ybar[pval > alpha]
	}
	options(warn = warn_state)
	xbar-ybar
}




##-----M.Padilla 2014

k.basic_BF_estimator <-function(x,y=NULL,...){
	if(length(y)>0){
		suppressWarnings(zA<-vapply(1:nrow(x),FUN.VALUE=numeric(1),FUN=function(i){
			z1<-try(logBF2posprob(BF.one.covariate(x=x[i,],y=y[i,],...)))
			if(is_err(z1)){NA} else {z1}
			}))
		}
	else{
		zA<-vapply(1:nrow(x),FUN.VALUE=numeric(1),FUN=function(i){
			z1<-try(logBF2posprob(BF.one.covariate(x=x[i,],...)))
			if(is_err(z1)){NA} else {z1}})
		}
	names(zA)<-rownames(x)
	zA
	}
k.AICc_estimator <-function(x,y=NULL,...){
	if(length(y)>0){
		suppressWarnings(zA<-vapply(1:nrow(x),FUN.VALUE=numeric(1),FUN=function(i){
			z1<-try(logBF2posprob(.5*AICc.one.covariate(x=x[i,],y=y[i,],...)))
			if(is_err(z1)){NA} else {z1}
			}))
		}
	else{
		zA<-vapply(1:nrow(x),FUN.VALUE=numeric(1),FUN=function(i){
			z1<-try(logBF2posprob(.5*AICc.one.covariate(x=x[i,],...)))
			if(is_err(z1)){NA} else {z1}})
		}
	names(zA)<-rownames(x)
	zA
	}
k.BIC_estimator <-function(x,y=NULL,...){
	if(length(y)>0){
		suppressWarnings(zA<-vapply(1:nrow(x),FUN.VALUE=numeric(1),FUN=function(i){
			z1<-try(logBF2posprob(.5*BIC.one.covariate(x=x[i,],y=y[i,],...)))
			if(is_err(z1)){NA} else {z1}
			}))
		}
	else{
		zA<-vapply(1:nrow(x),FUN.VALUE=numeric(1),FUN=function(i){
			z1<-try(logBF2posprob(.5*BIC.one.covariate(x=x[i,],...)))
			if(is_err(z1)){NA} else {z1}})
		}
	names(zA)<-rownames(x)
	zA
	}
	


k.nEA<- function(x,y=NULL,method=c("TWER","LFDR","FWER"),alpha=0.05,pval.fun=t.test_safe,arglis.pvalfun=list(),alternative,nulltype=0,...){
	paired<-F;err.num<-as.numeric(NA);param0<-0;param<-NULL;logx<-TRUE#if(method == "LFDR"){alternative <- "less"}
	
	z<-nprnSet2matrix(x,y=y,paired=paired)
	x<-z$x;y<-z$y;method<-method[1]
	param<-interest.param.fmeans(x=x,y=y,paired=paired,param=param,param0=param0,logx=logx)
	xbar<-param$xbar
	ybar<-param$ybar
	param<-param$param
	
	if(missing(alternative)&&method %in% c("TWER","FWER")){alternative<-"two.sided"}
	else if(missing(alternative)&&method %in% c("LFDR")){alternative<-"greater"}
	
	pval.x <- pval<-do.call(get_pvalues,c(list(x,y=y,test.fun=pval.fun,alternative=alternative),arglis.pvalfun))
	orig.pval<-pval
	pval<-pval[is.finite(pval)]
	if(length(pval)==0){zo<-param;return(zo)}
	err.pval<-rep(err.num,length(orig.pval))
	names(err.pval)<-names(orig.pval)
	
	if(method %in% c("FWER",'fwer')){
		z <- try(mt.rawp2adjp(pval, proc = "SidakSD",...))
		if(is_err(pval.x)){pval.x <- rep(err.num,length(pval));names(pval.x)<-names(stat)}
		else{
			pval.x<-z$adjp[,2]
			names(pval.x)<-names(pval)[z$index]}
		}
	else if(method %in% c("LFDR",'lfdr')){
		stat<- qnorm(pval);stat<-stat[is.finite(stat)]
		pval.x <- try(locfdr.safe(zz=stat, nulltype = nulltype,...))
		if(is_err(pval.x)){pval.x <- rep(err.num,length(stat))}
		names(pval.x)<-names(stat)}
	else if(method%in% c("TWER",'twer')){pval.x<-npval<-pval}
	else {stop('bad input method')}
	
	if(is_err(pval.x)||all(is.na(pval.x))){
		message(paste(method,'failed'))
		zo<-sameAsY(x=param,y=orig.pval);return(zo)}
	else{
		pval.x<-pval.x[is.finite(pval.x)]
		stopifnot(length(pval.x)>0)
		zo<-NULL}
	
	nxbar<-sameAsY(y=pval.x,x=xbar)

	if(length(y)==0){nxbar[pval.x > alpha] <- 0;nybar<-0}
	else{
		nybar<-sameAsY(y=pval.x,x=ybar)
		stopifnot(identical(names(nybar),names(nxbar)))
		nxbar[pval.x > alpha] <- nybar[pval.x > alpha]
		}

	sameAsY(y=orig.pval,x=(nxbar-nybar))

}
#


nscottberger.est<-function(x,y=NULL,scaled=F,...){
	z<-k.xy.prnset2matrix(x=x,y=y,paired=F)
	x<-z$x;y<-z$y
	err<-rep(as.numeric(NA),nrow(x))
	names(err)<-rownames(x)

	if(!scaled){zo<-try(scottbergermodel_estimator(x=x,y=y,...))}
	else if(scaled){zo<-try(scottbergermodel_scaled_data_estimator(x=x,y=y,...))}

	if(is_err(zo)){
		message(paste('scottberger estimator failed'))
		zo<-err}

	assert.is(zo,'numeric')
	zo<-MakeNames(zo,nmvar='X')
	sameAsY(x=zo,y=err)
}
#----
nIC.est<- function(x,y=NULL,opt=c('BF','AIC','BIC'),param0=NULL,param=NULL,logx=TRUE,...){
	opt<-opt[1];paired<-F;
	z<-nprnSet2matrix(x=x,y=y,paired=paired)
	x<-z$x;y<-z$y
	zp<-interest.param.fmeans(x=x,y=y,paired=paired,param=param,param0=param0,logx=logx)
	param<-zp$param
	param0<-zp$param0
	if(opt%in%c('basicBFestimator','BF','bf')){
		opt<-'basicBFestimator'
		p0<-k.basic_BF_estimator(x=x,y=y,...)
	}
	else if(opt%in%c('AICc.estimator','AIC','aic')){
		opt<-'AICc.estimator'
		p0<-k.AICc_estimator(x=x,y=y,...)
	}
	else if(opt%in%c('BICestimator','BIC','bic')){
		opt<-'BICestimator'
		p0<-k.BIC_estimator(x=x,y=y,...)
	}
	
	if(is_err(p0)){
		message(paste(opt,' failed'))
		p0<-rep(as.numeric(NA),length(param))
		names(p0)<-names(param)}
	
	p0<-sameAsY(x=p0,y=param)	
	p0*param0+(1 - p0)*param
}
nBF_estimator<-function(x,y=NULL,param0=NULL,param=NULL,logx=TRUE,...){
	nIC.est(x=x,y=y,opt='bf',param0=param0,param=param,logx=logx,...)
}
nAICc_estimator<-function(x,y=NULL,param0=NULL,param=NULL,logx=TRUE,...){
	nIC.est(x=x,y=y,opt='aic',param0=param0,param=param,logx=logx,...)
}
nBIC_estimator <-function(x,y=NULL,param0=NULL,param=NULL,logx=TRUE,...){
	nIC.est(x=x,y=y,opt='bic',param0=param0,param=param,logx=logx,...)
}
#
k.get.nulltype.from.opt<-function(opt){#opt<-"LFDR1" or "lfdr1", number in 0:3
	lfdrok<-grepl(pattern="lfdr", x=opt, ignore.case = TRUE,fixed=FALSE)
	zo<-list(lfdr.ok=lfdrok,nulltype=0,opt=opt)
	if(!lfdrok){return(zo)}
	zaux<-strsplit(x=opt,split=NULL)[[1]]
	suppressWarnings(nulltype<-as.numeric(zaux[length(zaux)]))
	if(!is.finite(nulltype)){nulltype<-0}
	stopifnot(is(nulltype,"numeric"))
	stopifnot(nulltype%in%c(0:3))
	opt<-"LFDR"
	zo$nulltype<-nulltype
	zo$opt<-"LFDR"
	zo
}
nhard.threshold.est<-function(x,y=NULL,opt=c('fold.change',"twer","fwer","lfdr", "lfdr0", "lfdr1"),alpha=0.05,
	pval.fun=t.test,arglis.pvalfun=list(),alternative="two.sided", ...){opt<-opt[1]
	zlfdrok<-k.get.nulltype.from.opt(opt);
	
	if(opt %in% c('fc','FC','fold.change')){
		z<-k.matrix2prnSet(x=x,y=y);x<-z$x;y<-z$y
		zo<-try(estimator_func_generator(probabilityNull = pseudoFdr, foldchange.threshold = alpha,...)(x=x,y=y))
	}
	else if(opt %in% c("TWER","FWER","twer","fwer","lfdr","LFDR") || zlfdrok$lfdr.ok){
		nulltype<-zlfdrok$nulltype
		opt<-zlfdrok$opt

		z<-k.xy.prnset2matrix(x=x,y=y,paired=F);x<-z$x;y<-z$y
		zo<-k.nEA(x=x,y=y,method=toupper(opt),alpha=alpha,pval.fun=pval.fun,
			arglis.pvalfun=arglis.pvalfun,alternative=alternative,nulltype=nulltype,...)
		}
	else(stop('bad opt'))
	err<-rep(as.numeric(NA),nrow(x))
	names(err)<-rownames(x)
	if(is_err(zo)){message(paste(opt,' failed'));zo<-err}

	assert.is(zo,'numeric')
	zo<-MakeNames(zo,nmvar='X')
	sameAsY(x=zo,y=err)
}
nfc.threshold.est<-function(x,y=NULL,threshold=0.5,...){nhard.threshold.est(x=x,y=y,opt=c('fold.change'),alpha=threshold,...)}
#

other.est<-function(x,y=NULL,opt=c("limma","pseudo","lfdr0","lfdr1"),pval.fun=t.test,alternative="greater",arglis.pvalfun=list(),...){opt<-opt[1]#,'locfdr.theonull','locfdr.empnull'
	
	moderated.t<-F;pfun<-NULL;param0<-0
	if(is(pval.fun,"character")){
		pval.fun<-pval.fun[1]
		if(pval.fun%in%c('mod.t','sFdr')){moderated.t<-T;pfun<-sFdr}}#probabilityNull <- eFdr0
	if(length(pfun)==0){
		pfun<-get_pvalfun(pval.fun=pval.fun,alternative=alternative,arglis.pvalfun=arglis.pvalfun)
		}
	z<-k.matrix2prnSet(x=x,y=y);x<-z$x;y<-z$y
	err<-rep(as.numeric(NA),nrow(x))
	names(err)<-rownames(x)
	
	if(opt%in%c('limma')){zo<-try(estimator_func_generator(probabilityNull = sFdr,...)(x=x,y=y))}
	else if(opt%in%c('pseudo')){zo<-try(estimator_func_generator(probabilityNull = pseudoFdr,...)(x=x,y=y))}
	else if(opt%in%c('lfdr0')){zo<-try(estimator_func_generator(probabilityNull = eFdr0, PValue = pfun,...)(x=x,y=y))}
	else if(opt%in%c('lfdr1')){zo<-try(estimator_func_generator(probabilityNull = eFdr, PValue = pfun,...)(x=x,y=y))}
	else{stop('bad opt')}
	
	if(is_err(zo)){message(paste(opt,' failed'));zo<-err}
	assert.is(zo,'numeric')
	zo<-MakeNames(zo,nmvar='X')
	sameAsY(x=zo,y=err)
}
nlimma.est<-function(x,y=NULL,...){
	other.est(x=x,y=y,opt=c('limma'),pval.fun=t.test,alternative='greater',arglis.pvalfun=list(),...)}
npseudo.est<-function(x,y=NULL,...){
	other.est(x=x,y=y,opt=c('pseudo'),pval.fun=t.test,alternative='greater',arglis.pvalfun=list(),...)}
nlfdr0.est<-function(x,y=NULL,pval.fun=t.test,alternative='greater',arglis.pvalfun=list(),...){
	other.est(x=x,y=y,opt=c('lfdr0'),pval.fun=pval.fun,alternative=alternative,arglis.pvalfun=arglis.pvalfun,...)}
nlfdr1.est<-function(x,y=NULL,pval.fun=t.test,alternative='greater',arglis.pvalfun=list(),...){
	other.est(x=x,y=y,opt=c('lfdr1'),pval.fun=pval.fun,alternative=alternative,arglis.pvalfun=arglis.pvalfun,...)}
#--------
k.locfdr.stat<- function(stat=NULL,pvalue=NULL,param=NULL,param0=NULL,nulltype=0, ...){
		
	z<-k.checkings.locfdr(pvalue=pvalue,stat=stat,param=param,nulltype=nulltype)
	param<-z$param
	stopifnot(is(param0,"numeric")&&length(param0)>0&&is.finite(param0))
	p0<-try(locfdr.safe(zz=z$stat,nulltype=z$nulltype, ...))

	if(is_err(p0)){p0<-rep(as.numeric(NA),length(param));names(p0)<-names(param)}
	else{if(length(names(p0))==0){names(p0)<-names(z$stat)}}
		 identical(names(z$stat),names(p0))
	p0<-sameAsY(x=p0,y=param)	
	zo<-p0*param0+(1 - p0)*param
	sameAsY(x=zo,y=z$orig.stat)

}


nlocfdr.est<-function(x.stat=NULL,y.pvalue=NULL,pval.fun=t.test,arglis.pvalfun=list(),alternative="greater",param=NULL,param0=NULL,logx=TRUE,nulltype=0,q.norm=T, ...){
	if(is(x.stat,"numeric")){stat<-x.stat;x<-NULL}
	else {x<-x.stat;stat<-NULL}
	if(is(y,"numeric")){pvalue<-y.pvalue;y<-NULL}
	else {y<-y.pvalue;pvalue<-NULL}

	stopifnot(length(x)>0||length(param)>0)
	stopifnot(length(x)>0||(length(stat)>0||length(pvalue)>0))
	
	if(length(x)>0){
		z<-k.xy.prnset2matrix(x=x,y=y,paired=F)
		x<-z$x;y<-z$y}
	
	if(length(param)==0||length(param0)==0){
		zp<-interest.param.fmeans(x=x,y=y,paired=F,param=param,param0=param0,logx=logx)
		if(length(param)==0){param<-zp$param}
		if(length(param0)==0){param0<-zp$param0}}
	
	if(length(stat)>0&&length(pvalue)>0&&q.norm){stat<-NULL}
	if(length(stat)>0&&length(pvalue)>0&&q.norm){pvalue<-NULL}
		
	if(length(stat)==0&&length(pvalue)==0&&q.norm){
		pvalue<-do.call(get_pvalues,c(list(x=x,y=y,test.fun=pval.fun,alternative=alternative),arglis.pvalfun))
		}
	else if(length(stat)==0&&length(pvalue)==0&&!q.norm){
		stat<-do.call(get_stats,c(list(x=x,y=y,test.fun=pval.fun),arglis.pvalfun))
		}
	if(length(stat)==0&&length(pvalue)>0&&q.norm){stat<-qnorm(pvalue)}
	else if (length(stat)==0&&length(pvalue)>0&&!q.norm){stop("error: q.norm=F and stat is empty")}		
				
	k.locfdr.stat(stat=stat,pvalue=pvalue,param=param,param0=param0,nulltype=nulltype, ...)
	
}
nlocfdr.x<-function(x=NULL,y=NULL,pval.fun=t.test,arglis.pvalfun=list(),alternative="greater",param=NULL,param0=NULL,logx=TRUE,nulltype=0,q.norm=T, ...){
	nlocfdr.est(x.stat=x,y.pvalue=y,pval.fun=pval.fun,arglis.pvalfun=arglis.pvalfun,alternative=alternative,param=param,param0=param0,nulltype=nulltype,q.norm=q.norm,logx=logx, ...)
}
nlocfdr.stat<-function(stat=NULL,pvalue=NULL,param=NULL,param0=0,nulltype=0,q.norm=T, ...){
	nlocfdr.est(x.stat=stat,y.pvalue=pvalue,param=param,param0=param0,nulltype=nulltype,q.norm=q.norm, ...)
}

#


