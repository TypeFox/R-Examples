dglars_ccd <- function(n,p,X,y,family,setting){
	storage.mode(n) <- "integer"
	storage.mode(p) <- "integer"
	storage.mode(X) <- "double"
	storage.mode(y) <- "double"
	g_hat <- as.double(2)
	b <- double((p+1)*setting$np)
	dev <- double(setting$np)
	A <- as.integer(1:p)
	nav <- as.integer(0)
	df <- integer(setting$np)
	g_seq <- double(setting$np)
	storage.mode(setting$np) <- "integer"
	storage.mode(setting$g0) <- "double"
	storage.mode(setting$nccd) <- "integer"
	storage.mode(setting$eps) <- "double"
	mthd <- ifelse(setting$method=="dgLASSO",1L,0L)
	conv <- integer(1)
	fit=switch(family,
			binomial=.Fortran("dglars_ccd_b",n=n,p=p,X=X,y=y,np=setting$np,g0=setting$g0,g_hat=g_hat,
							  nstp=setting$nccd,eps=setting$eps,mthd=mthd,b=b,dev=dev,g_seq=g_seq,
							  A=A,df=df,nav=nav,conv=conv),
			poisson=.Fortran("dglars_ccd_p",n=n,p=p,X=X,y=y,np=setting$np,g0=setting$g0,g_hat=g_hat,
							  nstp=setting$nccd,eps=setting$eps,mthd=mthd,b=b,dev=dev,g_seq=g_seq,
							  A=A,df=df,nav=nav,conv=conv)
			)
	fit <- make_dglars(fit,setting)
	fit
}