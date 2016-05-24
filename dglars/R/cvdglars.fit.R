cvdglars.fit <- function(X,y,family=c("binomial","poisson"),control=list()){
	this.call <- match.call()
	family <- match.arg(family)
	if(is.data.frame(X)) X <- as.matrix(X)
	if(is.null(colnames(X))) colnames(X) <- paste("x",1:dim(X)[2],sep="")
	Xdim <- dim(X)
	n <- Xdim[1]
	p <- Xdim[2]
	min_np <- min(n-1,p)
	if(length(y)!=n)
	stop("length(y)!=dim(X)[1]")
	if(is.character(y)) y <- factor(y)
	if(is.factor(y)){
		if(nlevels(y)!=2)
		stop("only factors with two levels are allowed")
		ny <- as.numeric(y) - 1
	} else ny <- y
	if(!is.list(control))
	stop("control is not a list")
	setting <- list(algorithm="pc",method="dgLASSO",nfold=10,foldid=NULL,ng=100,
					nv=min_np,np=NULL,g0=ifelse(p<n,1.0e-04,0.05),dg_max=0,nNR=50,
					NReps=1.0e-06,ncrct=50,cf=0.5,nccd=1.0e+05,eps=1.0e-05)
	nmsSetting <- names(setting)
	setting[(nms <- names(control))] <- control
	if(length(noNms <- nms[!nms %in% nmsSetting]))
	warning("unknow names in control: ",paste(noNms,collapse =", "))
	if(!setting$algorithm %in% c("pc","ccd"))
	stop("'algorithm' should be one of \"pc\", \"ccd\"")
	if(!setting$method %in% c("dgLASSO","dgLARS"))
	stop("'method' should be one of \"dgLASSO\",\"dgLARS\"")
	if(setting$nfold<1 | setting$nfold>n)
	stop("'nfold' should be an integer between 1 and n")
	if(is.null(setting$foldid)) setting$foldid=sample(n)
	else{
		if(min(setting$foldid)<1 | max(setting$foldid)>n)
		stop("Read the documentation for 'foldid' for more details")
	}
	if(setting$ng<1)
	stop("'ng' should be a non-negative integer")
	if(setting$nv<0 | setting$nv>min_np) 
	stop("'nv' should be an integer between 1 and min(n-1,p)")
	if(is.null(setting$np))
	setting$np <- ifelse(setting$algorithm=="pc",min_np*50,100L)
	if(setting$np<=1)
	stop("'np' should be a non-negative integer. Read the documentation more details")
	if(setting$g0<0)
	stop("'g0' should be a non-negative value. Read the documentation more details")
	if(setting$dg_max<0)
	stop("'dg_max' should be a non-negative value. Read the documentation more details")
	if(setting$eps<=0)
	stop("'eps' should be a non-negative value. Read the documentation more details")
	if(setting$ncrct<=0)
	stop("'ncrct' should be a non-negative value")
	if(setting$NReps<=0)
	stop("'NReps' should be a non-negative value")
	if(setting$nNR<=0)
	stop("'nNR' should be a non-negative integer")
	if(setting$cf<0 | setting$cf>1)
	stop("'cf' should be a value in the interval (0,1)")
	if(setting$nccd<=0)
	stop("'nccd' should be a non-negative integer")
	fit=switch(setting$algorithm,
			   pc=cvdglars_pc(n,p,X,ny,family,setting),
			   ccd=cvdglars_ccd(n,p,X,ny,family,setting)
			   )
	fit$call <- this.call
	fit$family <- family
	fit$y <- y
	if(fit$conv!=0) stop("dglars with algorithm ",fit$control$algorithm," does not converge with exit ",fit$conv)
	fit
}
