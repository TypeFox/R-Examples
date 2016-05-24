cvdglars_ccd <- function(n,p,X,y,family,setting){
	g <- (setting$ng-1):0/(setting$ng-1)
	b <- double(p+1)
	dev_m <- double(setting$ng)
	dev_v <- double(setting$ng)
	storage.mode(n) <- "integer"
	storage.mode(p) <- "integer"
	storage.mode(X) <- "double"
	storage.mode(y) <- "double"
	storage.mode(setting$foldid) <- "integer"
	storage.mode(setting$nfold) <- "integer"
	storage.mode(setting$ng) <- "integer"
	storage.mode(g) <- "double"
	g_hat <- double(1)
	mthd <- ifelse(setting$method=="dgLASSO",1L,0L)
	storage.mode(setting$g0) <- "double"
	storage.mode(setting$eps) <- "double"
	storage.mode(setting$np) <- "integer"
	storage.mode(setting$nccd) <- "integer"
	conv <- integer(1)
	fit=switch(family,
			binomial=.Fortran("cvdglars_ccd_b",n=n,p=p,X=X,y=y,foldid=setting$foldid,nfold=setting$nfold,
							  ng=setting$ng,g=g,b=b,dev_m=dev_m,dev_v=dev_v,g_hat=g_hat,mthd=mthd,g0=setting$g0,
							  eps=setting$eps,np=setting$np,nstp=setting$nccd,conv=conv),
			poisson=.Fortran("cvdglars_ccd_p",n=n,p=p,X=X,y=y,foldid=setting$foldid,nfold=setting$nfold,
							 ng=setting$ng,g=g,b=b,dev_m=dev_m,dev_v=dev_v,g_hat=g_hat,mthd=mthd,g0=setting$g0,
							 eps=setting$eps,np=setting$np,nstp=setting$nccd,conv=conv)
			)
	fit <- make_cvdglars(fit,setting)
	fit
}