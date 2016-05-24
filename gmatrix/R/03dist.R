# Distribution Functions
# 
# Author: nmorris
###############################################################################


###########################
#         normal
###########################
grnorm=function(n, mean = 0, sd = 1, type="d") {
	typeno=.type_num(type)
	if(typeno > 1L)
		stop("Normal variates must be of type 'double' or 'single.'")
	n=as.integer(n)[1]
	if(class(mean)!="gvector")
		mean=as.gvector(mean, type=typeno)
	if(class(sd)!="gvector")
		sd=as.gvector(sd, type=typeno)
	if(mean@type!=typeno)
		mean=convertType(mean,typeno)
	if(sd@type!=typeno)
		sd=convertType(sd,typeno)
	checkDevice(c(mean@device,sd@device))
	new("gvector",ptr=.Call("gpu_rnorm", n,mean@ptr,sd@ptr,mean@length,sd@length, typeno),length=n, type=typeno)
}



gdnorm = function (x, mean = 0, sd = 1, log = FALSE, type="d") {
	typeno=.type_num(type)
	if(typeno > 1L)
		stop("Normal variates must be of type 'double' or 'single.'")
	n=as.integer(length(x))
	
	if(class(mean)!="gvector")
		mean=as.gvector(mean, type=typeno)
	if(class(sd)!="gvector")
		sd=as.gvector(sd, type=typeno)
	if(class(x)!="gvector")
		x=as.gvector(x, type=typeno)
	
	if(mean@type!=typeno)
		mean=convertType(mean,typeno)
	if(sd@type!=typeno)
		sd=convertType(sd,typeno)
	if(x@type!=typeno)
		x=convertType(x,typeno)
	
	log=as.logical(log[1])
	if(log!=TRUE && log!=FALSE)
		stop("'log' must be TRUE or FALSE")
	
	#gpu_dnorm(SEXP in_n, SEXP in_x, SEXP in_mean, SEXP in_sd, SEXP in_n_mean, SEXP in_n_sd,
	#		SEXP in_log
	checkDevice(c(x@device, mean@device, sd@device))
	new("gvector",ptr=.Call("gpu_dnorm", n,x@ptr, mean@ptr,sd@ptr,mean@length,sd@length, log, typeno),length=n, type=typeno)
}


gqnorm = function (q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, warn=TRUE, type="d") {
	n=as.integer(length(q))
	
	typeno=.type_num(type)
	if(typeno > 1L)
		stop("Normal variates must be of type 'double' or 'single.'")
	
	if(class(mean)!="gvector")
		mean=as.gvector(mean, type=typeno)
	if(class(sd)!="gvector")
		sd=as.gvector(sd, type=typeno)
	if(class(q)!="gvector")
		q=as.gvector(q, type=typeno)
	
	if(mean@type!=typeno)
		mean=convertType(mean,typeno)
	if(sd@type!=typeno)
		sd=convertType(sd,typeno)
	if(q@type!=typeno)
		q=convertType(q,typeno)
	
	log.p=as.logical(log.p[1])
	if(log.p!=TRUE && log.p!=FALSE)
		stop("'log' must be TRUE or FALSE")
	lower.tail=lower.tail[1]
	if(lower.tail!=TRUE && lower.tail!=FALSE)
		stop("'log' must be TRUE or FALSE")
	if(warn)
		warning("'gqnorm' may not be accurate in the extreme tails.\n You can turn this message off with the 'warn' parameter.")
	#SEXP gpu_pnorm(SEXP in_n, SEXP in_x, SEXP in_mean, SEXP in_sd, SEXP in_n_mean, SEXP in_n_sd,
	#		SEXP in_log, SEXP in_lower);
	checkDevice(c(q@device, mean@device, sd@device))
	new("gvector",ptr=.Call("gpu_qnorm", n,q@ptr, mean@ptr,sd@ptr,mean@length,sd@length, log.p, lower.tail, typeno),length=n, type=typeno)
}


gpnorm = function (p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, warn=TRUE, type="d") {
	n=as.integer(length(p))
	typeno=.type_num(type)
	if(typeno > 1L)
		stop("Normal variates must be of type 'double' or 'single.'")
	
	if(class(mean)!="gvector")
		mean=as.gvector(mean, type=typeno)
	if(class(sd)!="gvector")
		sd=as.gvector(sd)
	if(class(p)!="gvector")
		p=as.gvector(p, type=typeno)
	
	if(mean@type!=typeno)
		mean=convertType(mean,typeno)
	if(sd@type!=typeno)
		sd=convertType(sd,typeno)
	if(p@type!=typeno)
		p=convertType(p,typeno)
	
	log.p=as.logical(log.p[1])
	if(log.p!=TRUE && log.p!=FALSE)
		stop("'log' must be TRUE or FALSE")
	lower.tail=lower.tail[1]
	if(lower.tail!=TRUE && lower.tail!=FALSE)
		stop("'log' must be TRUE or FALSE")
	if(warn)
		warning("'gpnorm' may not be accurate in the extreme tails.\n You can turn this message off with the 'warn' parameter.")
	#SEXP gpu_pnorm(SEXP in_n, SEXP in_x, SEXP in_mean, SEXP in_sd, SEXP in_n_mean, SEXP in_n_sd,
	#		SEXP in_log, SEXP in_lower);
	checkDevice(c(p@device, mean@device, sd@device))
	new("gvector",ptr=.Call("gpu_pnorm", n,p@ptr, mean@ptr,sd@ptr,mean@length,sd@length, log.p, lower.tail, typeno),length=n, type=typeno)
}
###########################
#         gamma
###########################
grgamma=function(n, shape, rate = 1, scale = 1/rate, type="d") {
	n=as.integer(n)[1]
	typeno=.type_num(type)
	if(typeno > 1L)
		stop("Normal variates must be of type 'double' or 'single.'")
	
	if(class(shape)!="gvector")
		shape=as.gvector(shape, type=typeno)
	if(class(scale)!="gmatrix")
		scale=as.gvector(scale, type=typeno)
	
	if(shape@type!=typeno)
		shape=convertType(shape,typeno)
	if(scale@type!=typeno)
		scale=convertType(scale,typeno)

	checkDevice(c(shape@device, scale@device))
	new("gvector",ptr=.Call("gpu_rgamma", n,shape@ptr,scale@ptr,shape@length,scale@length, typeno),length=n, type=typeno)
}

gdgamma=function(x, shape, rate = 1, scale = 1/rate, log = FALSE, type="d"){
	n=as.integer(length(x))
	typeno=.type_num(type)
	if(typeno > 1L)
		stop("Normal variates must be of type 'double' or 'single.'")
	
	if(class(shape)!="gvector")
		shape=as.gvector(shape, type=typeno)
	if(class(scale)!="gvector")
		scale=as.gvector(scale, type=typeno)
	if(class(x)!="gvector")
		x=as.gvector(x)
	
	if(shape@type!=typeno)
		shape=convertType(shape,typeno)
	if(scale@type!=typeno)
		scale=convertType(scale,typeno)
	if(x@type!=typeno)
		x=convertType(x,typeno)
	
	log=as.logical(log[1])
	if(log!=TRUE && log!=FALSE)
		stop("'log' must be TRUE or FALSE")
	
	checkDevice(c(x@device, shape@device, scale@device))
	new("gvector",ptr=.Call("gpu_dgamma", n,x@ptr, shape@ptr,scale@ptr,shape@length,scale@length, log, typeno),length=n, type=typeno)
}


###########################
#         uniform
###########################
grunif=function(n, min=0, max=1, type="d")
	
{
	n=as.integer(n)[1]
	typeno=.type_num(type)
	if(typeno>1L)
		stop("Ivalid type.")
	if(typeno > 1)
		stop("Normal variates must be of type 'double' or 'single.'")
	
	if(class(min)!="gvector")
		min=as.gvector(min, type=typeno)
	if(class(max)!="gmatrix")
		max=as.gvector(max, type=typeno)
	
	if(min@type!=typeno)
		min=convertType(min,typeno)
	if(min@type!=typeno)
		min=convertType(min,typeno)

	checkDevice(c(min@device, max@device))
	new("gvector",ptr=.Call("gpu_runif", n,min@ptr,max@ptr,min@length,max@length,typeno),
			length=n, type=typeno)
}

gdunif=function(x, min=0, max=1, log = FALSE, type="d"){
	n=as.integer(length(x))
	typeno=.type_num(type)
	if(typeno>1L)
		stop("Ivalid type.")
	if(typeno > 1)
		stop("Normal variates must be of type 'double' or 'single.'")
	
	if(class(min)!="gvector")
		min=as.gvector(min, type=typeno)
	if(class(max)!="gvector")
		max=as.gvector(max, type=typeno)
	if(class(x)!="gvector")
		x=as.gvector(x)
	
	if(min@type!=typeno)
		min=convertType(min,typeno)
	if(min@type!=typeno)
		min=convertType(min,typeno)
	if(x@type!=typeno)
		x=convertType(x,typeno)
	
	log=as.logical(log[1])
	if(log!=TRUE && log!=FALSE)
		stop("'log' must be TRUE or FALSE")
	checkDevice(c(x@device,min@device,max@device))
	new("gvector",ptr=.Call("gpu_dunif", n,x@ptr, min@ptr,max@ptr,min@length,max@length, log,typeno),length=n,type=typeno)
}


###########################
#         beta
###########################
grbeta=function(n, shape1, shape2, ncp = 0, type="d") {
	n=as.integer(n)[1]
	typeno=.type_num(type)
	if(typeno>1L)
		stop("Ivalid type.")
	if(typeno > 1)
		stop("Normal variates must be of type 'double' or 'single.'")
	
	if(!is.double(ncp))
		ncp=as.double(ncp)
	if(!all(ncp==0))
		stop("Noncental beta not implement on the GPU")
	
	if(class(shape1)!="gvector")
		shape1=as.gvector(shape1, type=typeno)
	if(class(shape2)!="gvector")
		shape2=as.gvector(shape2, type=typeno)

	if(shape1@type!=typeno)
		shape1=convertType(shape1,typeno)
	if(shape2@type!=typeno)
		shape2=convertType(shape2,typeno)

	checkDevice(c(shape1@device,shape2@device))
	new("gvector",ptr=.Call("gpu_rbeta", n,shape1@ptr,shape2@ptr,shape1@length,shape2@length,typeno),
			length=n, type=typeno)
}

gdbeta=function(x, shape1, shape2, ncp = 0, log = FALSE, type="d") {
	n=as.integer(length(x))
	typeno=.type_num(type)
	if(typeno>1L)
		stop("Ivalid type.")
	if(typeno > 1)
		stop("Normal variates must be of type 'double' or 'single.'")
	
	if(!all(ncp==0))
		stop("Noncental beta not implement on the GPU")
	if(class(shape1)!="gvector")
		shape1=as.gvector(shape1)
	if(class(shape2)!="gvector")
		shape2=as.gvector(shape2)
	if(class(x)!="gvector")
		x=as.gvector(x)
	
	if(shape1@type!=typeno)
		shape1=convertType(shape1,typeno)
	if(shape2@type!=typeno)
		shape2=convertType(shape2,typeno)
	if(x@type!=typeno)
		x=convertType(x,typeno)
	
	
	log=as.logical(log[1])
	if(log!=TRUE && log!=FALSE)
		stop("'log' must be TRUE or FALSE")
	checkDevice(c(x@device,shape1@device,shape2@device))
	new("gvector",ptr=.Call("gpu_dbeta", n,x@ptr, shape1@ptr,shape2@ptr,shape1@length,shape2@length, log,typeno),length=n,type=typeno)
}

###########################
#         binomial
###########################
grbinom=function(n, size, prob)
{
	n=as.integer(n)[1]
	if(class(size)!="gvector")
		size=as.gvector(size)
	if(class(prob)!="gmatrix")
		prob=as.gvector(prob)
	if(size@type==1L && prob@type==1L)
		typeno=1L
	else
		typeno=0L
	
	if(size@type!=typeno)
		size=convertType(size,typeno)
	if(prob@type!=typeno)
		prob=convertType(prob,typeno)
	
	checkDevice(c(size@device,prob@device))
	new("gvector",ptr=.Call("gpu_rbinom", n,size@ptr,prob@ptr,size@length,prob@length, typeno),length=n, type=2L)
}

gdbinom=function(x, size, prob, log = FALSE, type="d") {
	n=as.integer(length(x))
	typeno=.type_num(type)
	if(typeno>1L)
		stop("Ivalid type.")
	if(class(size)!="gvector")
		size=as.gvector(size)
	if(class(prob)!="gvector")
		prob=as.gvector(prob)
	if(class(x)!="gvector")
		x=as.gvector(x)
	
	log=log[1]
	if(log!=TRUE && log!=FALSE)
		stop("'log' must be TRUE or FALSE")
	
	if(size@type!=typeno)
		size=convertType(size,typeno)
	if(prob@type!=typeno)
		prob=convertType(prob,typeno)
	if(x@type!=typeno)
		x=convertType(x,typeno)
	
	checkDevice(c(x@device,size@device,prob@device))
	new("gvector",ptr=.Call("gpu_dbinom", n,x@ptr, size@ptr,prob@ptr,size@length,prob@length, log, typeno),length=n, type=typeno)
}
###########################
#         poison
###########################
grpois=function(n, lambda)  {
	n=as.integer(n)[1]
	if(class(lambda)!="gvector")
		lambda=as.gvector(lambda)
	if(lambda@type>1L)
		type(lambda)=0L
	checkDevice(c(lambda@device))
	new("gvector",ptr=.Call("gpu_rpois", n,lambda@ptr,lambda@length,lambda@type),
			length=n, type=2L)
}

gdpois=function(x, lambda, log = FALSE, type="d") {
	n=as.integer(length(x))
	typeno=.type_num(type)
	if(typeno>1L)
		stop("Ivalid type.")
	if(class(lambda)!="gvector")
		lambda=as.gvector(lambda)
	if(class(x)!="gvector")
		x=as.gvector(x)
	log=log[1]
	if(log!=TRUE && log!=FALSE)
		stop("'log' must be TRUE or FALSE")
	if(lambda@type!=typeno)
		lambda=convertType(lambda,typeno)
	if(x@type!=typeno)
		x=convertType(x,typeno)
	checkDevice(c(x@device,lambda@device))
	new("gvector",ptr=.Call("gpu_dpois", n, x@ptr, lambda@ptr,lambda@length, log, typeno),length=n, type=typeno)
}

###########################
#         rsample
###########################
rsample = function(P, log=TRUE) {
	if(class(P)!="gmatrix")
		stop("Object must be of class 'gmatrix.'")
	if(!log)
		P=log(P)
	#SEXP gpu_rsample(SEXP in_P, SEXP in_rows, SEXP in_cols, SEXP in_norm, SEXP in_type);
	norm = gRowLogSums(P)
	return(new("gvector", ptr=.Call("gpu_rsample",P@ptr, nrow(P), ncol(P),norm@ptr, P@type), length=nrow(P), type=2L))
}

