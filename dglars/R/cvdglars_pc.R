cvdglars_pc <- function(n,p,X,y,family,setting){
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
	storage.mode(setting$nv) <- "integer"
	mthd <- ifelse(setting$method=="dgLASSO",1L,0L)
	storage.mode(setting$g0) <- "double"
	storage.mode(setting$dg_max) <- "double"
	storage.mode(setting$eps) <- "double"
	storage.mode(setting$np) <- "integer"
	storage.mode(setting$ncrct) <- "integer"
	storage.mode(setting$cf) <- "double"
	storage.mode(setting$NReps) <- "double"
	storage.mode(setting$nNR) <- "integer"
	conv <- integer(1)
	fit=switch(family,
			binomial=.Fortran("cvdglars_pc_b",n=n,p=p,X=X,y=y,foldid=setting$foldid,nfold=setting$nfold,
							  ng=setting$ng,g=g,b=b,dev_m=dev_m,dev_v=dev_v,g_hat=g_hat,nv=setting$nv,mthd=mthd,
							  g0=setting$g0,dg_max=setting$dg_max,eps=setting$eps,np=setting$np,ncrct=setting$ncrct,
							  cf=setting$cf,NReps=setting$NReps,nNR=setting$nNR,conv=conv),
			poisson=.Fortran("cvdglars_pc_p",n=n,p=p,X=X,y=y,foldid=setting$foldid,nfold=setting$nfold,
							 ng=setting$ng,g=g,b=b,dev_m=dev_m,dev_v=dev_v,g_hat=g_hat,nv=setting$nv,mthd=mthd,
							 g0=setting$g0,dg_max=setting$dg_max,eps=setting$eps,np=setting$np,ncrct=setting$ncrct,
							 cf=setting$cf,NReps=setting$NReps,nNR=setting$nNR,conv=conv)
			   )
	fit <- make_cvdglars(fit,setting)
	fit
}