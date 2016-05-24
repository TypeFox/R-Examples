dglars_pc <- function(n,p,X,y,family,setting){
	storage.mode(n) <- "integer"
	storage.mode(p) <- "integer"
	storage.mode(X) <- "double"
	storage.mode(y) <- "double"
	b <- double((p+1)*setting$np)
	ru <- double(p*setting$np)
	dev <- double(setting$np)
	A <- as.integer(1:p)
	storage.mode(setting$nv) <- "integer"
	nav <- as.integer(0)
	df <- integer(setting$np)
	g_seq <- double(setting$np)	
	storage.mode(setting$g0) <- "double"
	g_hat <- as.double(2)
	storage.mode(setting$dg_max) <- "double"
	storage.mode(setting$eps) <- "double"
	storage.mode(setting$np) <- "integer"
	storage.mode(setting$ncrct) <- "integer"
	storage.mode(setting$cf) <- "double"
	storage.mode(setting$NReps) <- "double"
	storage.mode(setting$nNR) <- "integer"
	mthd <- ifelse(setting$method=="dgLASSO",1L,0L)
	conv <- integer(1)
	fit=switch(family,
			binomial=.Fortran("dglars_pc_b",n=n,p=p,X=X,y=y,b=b,ru=ru,dev=dev,A=A,
							  nv=setting$nv,nav=nav,df=df,g_seq=g_seq,mthd=mthd,
							  g0=setting$g0,g_hat=g_hat,dg_max=setting$dg_max,
							  eps=setting$eps,np=setting$np,ncrct=setting$ncrct,
							  cf=setting$cf,NReps=setting$NReps,nNR=setting$nNR,conv=conv),
			poisson=.Fortran("dglars_pc_p",n=n,p=p,X=X,y=y,b=b,ru=ru,dev=dev,A=A,
							 nv=setting$nv,nav=nav,df=df,g_seq=g_seq,mthd=mthd,
							 g0=setting$g0,g_hat=g_hat,dg_max=setting$dg_max,
							 eps=setting$eps,np=setting$np,ncrct=setting$ncrct,
							 cf=setting$cf,NReps=setting$NReps,nNR=setting$nNR,conv=conv),
			   )
	fit <- make_dglars(fit,setting)
	fit
}