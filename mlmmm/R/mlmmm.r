############################################################################
## R code: "mlmmm.r" which defines "mlmmm.em" and "mlmmm.em" #
############################################################################

###############################################################################
# MLE in multivariate linear mixed models with possible(!) missing values
###############################################################################
#####
#
mlmmm.em<-function(y,subj,pred,xcol,zcol,start,maxits=200,eps=.0001){
	if(any(is.na(pred)))
		stop("missing values in pred not allowed")
        # change y and pred to matrices, if necessary
	if(is.vector(y)) y<-matrix(y,ncol=1)
	if(is.vector(pred)) pred<-matrix(pred,ncol=1)
	m<-as.integer(length(table(subj)))
	ntot<-as.integer(nrow(y))
	
	nmax<-as.integer(max(table(subj)))
	r<-as.integer(ncol(y))
	p<-length(xcol)
	q<-length(zcol)
	ggs<-as.integer(round(((q*r)*((q*r)+1)/2)+r*(r+1)/2))
	zcol<-as.integer(zcol)
	xcol<-as.integer(xcol)
	pcol<-as.integer(ncol(pred))
	#
	{if(missing(start)){
		beta<-matrix(0,p,r)
		sigma<-matrix(0,r,r)
		psi<-matrix(0,q*r,q*r)
		epsi<-matrix(0,ntot,r)
		sflag<-as.integer(0)}
	 else{
		beta<-start$beta
		sigma<-start$sigma
		psi<-start$psi
		epsi<-matrix(0,ntot,r)
		sflag<-as.integer(1)
		storage.mode(eps)<-"double"
		storage.mode(beta)<-"double"
		storage.mode(sigma)<-"double"
		storage.mode(psi)<-"double"}}
	cat(" ###  performing mle for mlmm with NA values ### ")
	now<-proc.time()
	#
	# create rmat, npatt and patt to keep track of missingness patterns
	rmat<-1-1*is.na(y)
	storage.mode(rmat)<-"integer"
	revcpatt<-rep("",ntot)
        for(i in 1:r) revcpatt<-paste(as.character(rmat[,i]),revcpatt,sep="")
	nulpat0<-""
        nulpat2<-""
        for(i in 1:r){
           nulpat0<-paste(nulpat0,"0",sep="")
           nulpat2<-paste(nulpat2,"2",sep="")}
        revcpatt[revcpatt==nulpat0]<-nulpat2
	tmp<-rev(table(revcpatt))
	npatt<-length(tmp)
	if(any(revcpatt==nulpat2)) npatt<-npatt-1
	ww<-!duplicated(revcpatt)
	upatt<-revcpatt[ww]
	rmat<-rmat[ww,]
        if(r==1) rmat<-matrix(rmat,ncol=1)
	ww<-rev(order(upatt))
	upatt<-upatt[ww]
	rmat<-matrix(rmat,ncol=r,nrow=length(rev(order(upatt))))
	rmat<-rmat[ww,]
        if(r==1) rmat<-matrix(rmat,ncol=1)
	if(any(upatt==nulpat2)){
		rmat<-rmat[-1,]
		upatt<-upatt[-1]}
        patt<-integer(ntot)
	patt[revcpatt==nulpat2]<-0
	for(i in 1:npatt) patt[revcpatt==upatt[i]]<-i
	storage.mode(npatt)<-"integer"
	storage.mode(rmat)<-"integer"
	storage.mode(patt)<-"integer"
	iposn<-as.integer(1:ntot)
	ww<-order(patt)
	iposn<-iposn[ww]
	pstfin<-matrix(0,npatt,2)
	{if(any(patt==0)){
		sst<-tmp[1]+1
		for(i in 2:(npatt+1)){
			pstfin[i-1,1]<-sst
			pstfin[i-1,2]<-sst+tmp[i]-1
			sst<-sst+tmp[i]}}
	else{
		sst<-1
		for(i in 1:npatt){
			pstfin[i,1]<-sst
			pstfin[i,2]<-sst+tmp[i]-1
			sst<-sst+tmp[i]}}}
	storage.mode(pstfin)<-"integer"
	#
	storage.mode(y)<-"double"
	y[is.na(y)]<--999.99
	storage.mode(pred)<-"double"
	#####
	#####
	tmp<-.Fortran("mlmmem2",
		     intinput= as.integer(c(ntot,
		       m,
		       r,
		       p,
		       q,
		       subj,
		       nmax,
		       iposn,
		       npatt,
		       pstfin,
		       patt,
		       rmat,
		       pcol,
		       xcol,
		       zcol,
		       maxits,
		       ggs,
		       sflag)),
		     intoutpt= integer(4+3*m),
		     dbinput= as.double(c(pred,
		       y,
		       sigma,
		       beta,
		       psi,
		       eps,
		       epsi)),
		     dboutput= numeric(r*nmax*r*nmax*10),
	#
		     w=array(0,c(r*nmax,r*nmax,m)),

		     wkqb2=matrix(0,nmax,r),
		     vdel=numeric(r*nmax),
        #
		     uszxb=numeric(r*q),usotzo=matrix(0,r*q,r*nmax),
		     usotzm=matrix(0,r*q,r*nmax),wxbw=numeric(r*nmax),
		     wxbwo=numeric(r*nmax),wxbwm=numeric(r*nmax),
		     wkeb2=matrix(0,r*q,r*nmax),eb=matrix(0,r*q,m),
		     wxbeta=matrix(0,ntot,r),wxbetazeb=matrix(0,ntot,r),
		     varb=array(0,c(r*q,r*q,m)),wkrrpt=array(0,c(r,r,npatt)),
		     wkrrb21=array(0,c(r,r,npatt)),
	#
		     eystar=matrix(0,ntot,r),ey=matrix(0,ntot,r),
		     u=array(0,c(r*q,r*q,m)),
		     ztz=array(0,c(q,q,m)),
		     xtw=matrix(0,p*r,nmax*r),xtwx=matrix(0,p*r,p*r),
		     xtwy=numeric(p*r),xtwxinv=matrix(0,p*r,p*r),
		     wkqq1=matrix(0,r*q,r*q),wkqq2=matrix(0,r*q,r*q),
        #
		     wkqq3=matrix(0,r*q,r*q),wkrr1=matrix(0,r,r),
		     wkrr2=matrix(0,r,r),wksigtz=array(0,c(r*q,r*nmax,m)),
		     wkqqu=array(0,c(r*q,r*q,m)),
		     wkqnm=array(0,c(r*q,r*nmax,m)),
		     obeta=matrix(0,p,r),
		     osigma=matrix(0,r,r),opsi=array(0,c(r*q,r*q)),
		     llvec=numeric(as.integer(maxits)),
		     llovec=numeric(as.integer(maxits)),
		     wkg=rep(0,ggs),wkgg=matrix(0,ggs,ggs),wkpr=matrix(0,p,r),
		     wkpp=matrix(0,p,p),xtxinv=matrix(0,p,p))
        ######## contents of "intinput" #######
	in1 <- 1
	in2 <- in1 + 1
	in3 <- in2 + 1
	in4 <- in3 + 1
	in5 <- in4 + 1
	ntot <- tmp$intinput[in1]
	m <- tmp$intinput[in2]
	r <- tmp$intinput[in3]
	p <- tmp$intinput[in4]
	q <- tmp$intinput[in5]
	in6 <- in5 + 1
	in7 <- in6 + ntot
	in8 <- in7 + 1
	in9 <- in8 + ntot
	in10 <- in9 + 1
	in11 <- in10 + 2*npatt
	in12 <- in11 + ntot
	in13 <- in12 + r*npatt
	in14 <- in13 + 1
	in15 <- in14 + p
	in16 <- in15 + q
	in17 <- in16 + 1
	in18 <- in17 + 1
	subj <- tmp$intinput[in6:(in7-1)]
	nmax <- tmp$intinput[in7]
	iposn <- tmp$intinput[in8:(in9-1)]
	npatt <- tmp$intinput[in9]
	pstfin <- matrix(tmp$intinput[in10:(in11-1)],nrow=npatt)
	patt <- tmp$intinput[in11:(in12-1)]
	rmat <- matrix(tmp$intinput[in12:(in13-1)],nrow=npatt)
	pcol <- tmp$intinput[in13]
	xcol <- tmp$intinput[in14:(in15-1)]
	zcol <- tmp$intinput[in15:(in16-1)]
	maxits <- tmp$intinput[in16]
	ggs <- tmp$intinput[in17]
	sflag <- tmp$intinput[in18]
        ######## contents of "dbinput" #######
	isub0 <- r*q
	idi1 <- 1
	idi2 <- idi1 + ntot*pcol
	idi3 <- idi2 + r*ntot
	idi4 <- idi3 + r*r
	idi5 <- idi4 + p*r
	idi6 <- idi5 + isub0*isub0
	idi7 <- idi6 + 1
	pred <- matrix(tmp$dbinput[idi1:(idi2-1)],nrow= ntot)
	y <- matrix(tmp$dbinput[idi2:(idi3-1)],nrow= ntot)
	sigma <- matrix(tmp$dbinput[idi3:(idi4-1)],nrow=r)
	beta <- matrix(tmp$dbinput[idi4:(idi5-1)],nrow= p)
	psi <- array(tmp$dbinput[idi5:(idi6-1)],dim=c(isub0,isub0))
	eps <- tmp$dbinput[idi6]
	epsi <- matrix(tmp$dbinput[idi7:(idi7+ntot*r-1)],nrow= ntot)
        ######## contents of "intoutpt" ######
	io1 <- 1
	io2 <- io1 + m
	io3 <- io2 + m
	io4 <- io3 + 1
	io5 <- io4 + m
	io6 <- io5 + 1
	io7 <- io6 + 1
	ist <- tmp$intoutpt[io1]
	ifin <- tmp$intoutpt[io2:(io3-1)]
	nstar <- tmp$intoutpt[io3]
	nstari <- tmp$intoutpt[io4:(io5-1)]
	iter <- tmp$intoutpt[io5]
	msg <- tmp$intoutpt[io6]
	cvgd <- tmp$intoutpt[io7]
        ######## contents of "dboutput" ######
	isub1 <- r*nmax
	isub <- isub1*isub1
	ido1 <- 1
	ido2 <- ido1 + isub
	ido3 <- ido2 + isub
	ido4 <- ido3 + isub
	ido5 <- ido4 + isub
	ido6 <- ido5 + isub
	ido7 <- ido6 + isub
	ido8 <- ido7 + isub
	ido9 <- ido8 + isub
	ido10 <- ido9 + isub
	wo <-       matrix(tmp$dboutput[ido1:(ido2-1)],nrow= isub1)
	wo1 <-      matrix(tmp$dboutput[ido2:(ido3-1)],nrow= isub1)
	wm <-       matrix(tmp$dboutput[ido3:(ido4-1)],nrow= isub1)
	wom <-      matrix(tmp$dboutput[ido4:(ido5-1)],nrow= isub1)
	wkwmm1 <-   matrix(tmp$dboutput[ido5:(ido6-1)],nrow= isub1)
	wkwmm2 <-   matrix(tmp$dboutput[ido6:(ido7-1)],nrow= isub1)
	eyyt <-     matrix(tmp$dboutput[ido7:(ido8-1)],nrow= isub1)
	eyxyxt <-   matrix(tmp$dboutput[ido8:(ido9-1)],nrow= isub1)
	wkeyxyxt <- matrix(tmp$dboutput[ido9:(ido10-1)],nrow= isub1)
	wkqnm1 <-   matrix(tmp$dboutput[ido10:(ido10+isub-1)],nrow= isub1)
	clock<-proc.time()-now
	cat("\n")
	{if(msg==1)
	 warning("xtx is not full rank, failed for calculating beta<-(0)")
	else if(msg==2)
	 warning("Value of psi or sigma or U<-i
	          became non-pos.def.during iterations")
	else if(msg==3)
	 warning("Value of Var(y<-i(obs)) became non-pos.def.during iterations")
	else if(msg==4)
          warning("GLS failed for start vals, xtwx not full rank")
	else if(msg==5)
	  warning("Value of psi became non-pos.def. during iterations")
	else if(msg==6)
	  warning("Value of sigma became non-pos.def. during iterations")
        else if(msg==7)
          warning("log-density not concave at one or more scoring steps")}
        llvec<-tmp$llvec[1:iter]
	llovec<-tmp$llovec[1:iter]
	converged<-cvgd==as.integer(1)
	if(!converged) warning(paste("did not converge by",
	   format(iter),"iterations"))
	#
	list(beta=beta,sigma=sigma,psi=psi,eb=tmp$eb,varb=tmp$varb,xtwxinv=tmp$xtwxinv,
			converged=converged,iter=iter,npatt=npatt,pstfin=pstfin,iposn=iposn,patt=patt,rmat=rmat,
			logll=llvec,logoll=llovec,clock=clock)}
###########################################################################
###########################################################################
# Version of em for mlmm that allows psi to be block-diagonal.
# This function looks the same as the mlmem, except that psi is
# now an array of dim. (q x q x r). Also the workspaces wkqq1 and wkqq2
# are changed to array(0,c(q,q,r)) and matrix(0,q,q). Two additional
# workspaces are added: wkrqrq1(r*q,r*q) and wkrqrq2(r*q,r*q).
#
mlmmmbd.em<-function(y,subj,pred,xcol,zcol,start,maxits=100,eps=.01){
	if(any(is.na(pred)))
		stop("missing values in pred not allowed")
        # change y and pred to matrices, if necessary
	if(is.vector(y)) y<-matrix(y,ncol=1)
	if(is.vector(pred)) pred<-matrix(pred,ncol=1)
	m<-as.integer(length(table(subj)))
	ntot<-as.integer(nrow(y))
	##for now
	nmax<-as.integer(max(table(subj)))
	r<-as.integer(ncol(y))
	p<-length(xcol)
	q<-length(zcol)
	ggs<-as.integer(round((r*(q*(q+1))/2)+r*(r+1)/2))
	zcol<-as.integer(zcol)
	xcol<-as.integer(xcol)
	pcol<-as.integer(ncol(pred))
	#
	# for now starting values are assumed to be supplied
	{if(missing(start)){
		beta<-matrix(0,p,r)
		sigma<-matrix(0,r,r)
		psi<-array(0,c(q,q,r))
		epsi<-matrix(0,ntot,r)
		sflag<-as.integer(0)}
	 else{
		beta<-start$beta
		sigma<-start$sigma
		psi<-start$psi
		epsi<-matrix(0,ntot,r)
		sflag<-as.integer(1)
		storage.mode(psi)<-"double"
		storage.mode(beta)<-"double"
		storage.mode(sigma)<-"double"}}
	cat("performing block-diagonal version of EM in mlmm with NA values")
	now<-proc.time()
	#
	# create rmat, npatt and patt to keep track of missingness patterns
	rmat<-1-1*is.na(y)
	storage.mode(rmat)<-"integer"
	revcpatt<-rep("",ntot)
        for(i in 1:r) revcpatt<-paste(as.character(rmat[,i]),revcpatt,sep="")
	nulpat0<-""
        nulpat2<-""
        for(i in 1:r){
           nulpat0<-paste(nulpat0,"0",sep="")
           nulpat2<-paste(nulpat2,"2",sep="")}
        revcpatt[revcpatt==nulpat0]<-nulpat2
	tmp<-rev(table(revcpatt))
	npatt<-length(tmp)
	if(any(revcpatt==nulpat2)) npatt<-npatt-1
	ww<-!duplicated(revcpatt)
	upatt<-revcpatt[ww]
	rmat<-rmat[ww,]
        if(r==1) rmat<-matrix(rmat,ncol=1)
	ww<-rev(order(upatt))
	upatt<-upatt[ww]
	rmat<-matrix(rmat,ncol=r,nrow=length(rev(order(upatt))))
	rmat<-rmat[ww,]
        if(r==1) rmat<-matrix(rmat,ncol=1)
	if(any(upatt==nulpat2)){
		rmat<-rmat[-1,]
		upatt<-upatt[-1]}
        patt<-integer(ntot)
	patt[revcpatt==nulpat2]<-0
	for(i in 1:npatt) patt[revcpatt==upatt[i]]<-i
	storage.mode(npatt)<-"integer"
	storage.mode(rmat)<-"integer"
	storage.mode(patt)<-"integer"
	iposn<-as.integer(1:ntot)
	ww<-order(patt)
	iposn<-iposn[ww]
	pstfin<-matrix(0,npatt,2)
	{if(any(patt==0)){
		sst<-tmp[1]+1
		for(i in 2:(npatt+1)){
			pstfin[i-1,1]<-sst
			pstfin[i-1,2]<-sst+tmp[i]-1
			sst<-sst+tmp[i]}}
	else{
		sst<-1
		for(i in 1:npatt){
			pstfin[i,1]<-sst
			pstfin[i,2]<-sst+tmp[i]-1
			sst<-sst+tmp[i]}}}
	storage.mode(pstfin)<-"integer"
	#
	storage.mode(y)<-"double"
	y[is.na(y)]<--999.99
	storage.mode(pred)<-"double"
	#####
	#####
	tmp<-.Fortran("mlmmembd2",
		     intinput= as.integer(c(ntot,
		       m,
		       r,
		       p,
		       q,
		       subj,
		       nmax,
		       iposn,
		       npatt,
		       pstfin,
		       patt,
		       rmat,
		       pcol,
		       xcol,
		       zcol,
		       maxits,
		       ggs,
		       sflag)),
		     intoutpt= integer(4+3*m),
		     dbinput= as.double(c(pred,
		       y,
		       sigma,
		       beta,
		       psi,
		       eps,
		       epsi)),
		     dboutput= numeric(r*nmax*r*nmax*10),
	#
		     w=array(0,c(r*nmax,r*nmax,m)),
		     wkqb2=matrix(0,nmax,r),
		     vdel=numeric(r*nmax),
        #
		     uszxb=numeric(r*q),usotzo=matrix(0,r*q,r*nmax),
		     usotzm=matrix(0,r*q,r*nmax),wxbw=numeric(r*nmax),
		     wxbwo=numeric(r*nmax),wxbwm=numeric(r*nmax),
		     wkeb2=matrix(0,r*q,r*nmax),eb=matrix(0,r*q,m),
		     wxbeta=matrix(0,ntot,r),wxbetazeb=matrix(0,ntot,r),
		     varb=array(0,c(r*q,r*q,m)),wkrrpt=array(0,c(r,r,npatt)),
		     wkrrb21=array(0,c(r,r,npatt)),
	#
		     eystar=matrix(0,ntot,r),ey=matrix(0,ntot,r),
		     u=array(0,c(r*q,r*q,m)),
		     ztz=array(0,c(q,q,m)),
		     xtw=matrix(0,p*r,nmax*r),xtwx=matrix(0,p*r,p*r),
		     xtwy=numeric(p*r),xtwxinv=matrix(0,p*r,p*r),
		     wkrqrq1=matrix(0,r*q,r*q),wkrqrq2=matrix(0,r*q,r*q),
        #
		     wkqq1bd=array(0,c(q,q,r)),wkqq2bd=matrix(0,q,q),
		     wkqq3=matrix(0,r*q,r*q),wkrr1=matrix(0,r,r),
		     wkrr2=matrix(0,r,r),wksigtz=array(0,c(r*q,r*nmax,m)),
		     wkqqu=array(0,c(r*q,r*q,m)),
		     wkqnm=array(0,c(r*q,r*nmax,m)),
		     obeta=matrix(0,p,r),
		     osigma=matrix(0,r,r),opsi=array(0,c(q,q,r)),
                     #####don't know how to calculate llvec
		     llvec=numeric(as.integer(maxits)),
		     llovec=numeric(as.integer(maxits)),
		     wkg=rep(0,ggs),wkgg=matrix(0,ggs,ggs),wkpr=matrix(0,p,r),
		     wkpp=matrix(0,p,p),xtxinv=matrix(0,p,p))
######## contents of "intinput" #######
	in1 <- 1
	in2 <- in1 + 1
	in3 <- in2 + 1
	in4 <- in3 + 1
	in5 <- in4 + 1
	ntot <- tmp$intinput[in1]
	m <- tmp$intinput[in2]
	r <- tmp$intinput[in3]
	p <- tmp$intinput[in4]
	q <- tmp$intinput[in5]
	in6 <- in5 + 1
	in7 <- in6 + ntot
	in8 <- in7 + 1
	in9 <- in8 + ntot
	in10 <- in9 + 1
	in11 <- in10 + 2*npatt
	in12 <- in11 + ntot
	in13 <- in12 + r*npatt
	in14 <- in13 + 1
	in15 <- in14 + p
	in16 <- in15 + q
	in17 <- in16 + 1
	in18 <- in17 + 1
	subj <- tmp$intinput[in6:(in7-1)]
	nmax <- tmp$intinput[in7]
	iposn <- tmp$intinput[in8:(in9-1)]
	npatt <- tmp$intinput[in9]
	pstfin <- matrix(tmp$intinput[in10:(in11-1)],nrow=npatt)
	patt <- tmp$intinput[in11:(in12-1)]
	rmat <- matrix(tmp$intinput[in12:(in13-1)],nrow=npatt)
	pcol <- tmp$intinput[in13]
	xcol <- tmp$intinput[in14:(in15-1)]
	zcol <- tmp$intinput[in15:(in16-1)]
	maxits <- tmp$intinput[in16]
	ggs <- tmp$intinput[in17]
	sflag <- tmp$intinput[in18]
        ######## contents of "dbinput" #######
	idi1 <- 1
	idi2 <- idi1 + ntot*pcol
	idi3 <- idi2 + r*ntot
	idi4 <- idi3 + r*r
	idi5 <- idi4 + p*r
	idi6 <- idi5 + q*q*r
	idi7 <- idi6 + 1
	pred <- matrix(tmp$dbinput[idi1:(idi2-1)],nrow= ntot)
	y <- matrix(tmp$dbinput[idi2:(idi3-1)],nrow= ntot)
	sigma <- matrix(tmp$dbinput[idi3:(idi4-1)],nrow=r)
	beta <- matrix(tmp$dbinput[idi4:(idi5-1)],nrow= p)
	psi <- array(tmp$dbinput[idi5:(idi6-1)],dim=c(q,q,r))
	eps <- tmp$dbinput[idi6]
	epsi <- matrix(tmp$dbinput[idi7:(idi7+ntot*r-1)],nrow= ntot)
        ######## contents of "intoutpt" ######
	io1 <- 1
	io2 <- io1 + m
	io3 <- io2 + m
	io4 <- io3 + 1
	io5 <- io4 + m
	io6 <- io5 + 1
	io7 <- io6 + 1
	ist <- tmp$intoutpt[io1]
	ifin <- tmp$intoutpt[io2:(io3-1)]
	nstar <- tmp$intoutpt[io3]
	nstari <- tmp$intoutpt[io4:(io5-1)]
	iter <- tmp$intoutpt[io5]
	msg <- tmp$intoutpt[io6]
	cvgd <- tmp$intoutpt[io7]
        ######## contents of "dboutput" ######
	isub1 <- r*nmax
	isub <- isub1*isub1
	ido1 <- 1
	ido2 <- ido1 + isub
	ido3 <- ido2 + isub
	ido4 <- ido3 + isub
	ido5 <- ido4 + isub
	ido6 <- ido5 + isub
	ido7 <- ido6 + isub
	ido8 <- ido7 + isub
	ido9 <- ido8 + isub
	ido10 <- ido9 + isub
	wo <-       matrix(tmp$dboutput[ido1:(ido2-1)],nrow= isub1)
	wo1 <-      matrix(tmp$dboutput[ido2:(ido3-1)],nrow= isub1)
	wm <-       matrix(tmp$dboutput[ido3:(ido4-1)],nrow= isub1)
	wom <-      matrix(tmp$dboutput[ido4:(ido5-1)],nrow= isub1)
	wkwmm1 <-   matrix(tmp$dboutput[ido5:(ido6-1)],nrow= isub1)
	wkwmm2 <-   matrix(tmp$dboutput[ido6:(ido7-1)],nrow= isub1)
	eyyt <-     matrix(tmp$dboutput[ido7:(ido8-1)],nrow= isub1)
	eyxyxt <-   matrix(tmp$dboutput[ido8:(ido9-1)],nrow= isub1)
	wkeyxyxt <- matrix(tmp$dboutput[ido9:(ido10-1)],nrow= isub1)
	wkqnm1 <-   matrix(tmp$dboutput[ido10:(ido10+isub-1)],nrow= isub1)
	clock<-proc.time()-now
	cat("\n")
	{if(msg==1)
	 warning("xtx is not full rank, failed for calculating beta<-(0)")
	else if(msg==2)
	 warning("Value of psi or sigma or U<-i
	          became non-pos.def.during iterations")
	else if(msg==3)
	 warning("Value of Var(y<-i(obs)) became non-pos.def.during iterations")
	else if(msg==4)
          warning("GLS failed for start vals, xtwx not full rank")
	else if(msg==5)
	  warning("Value of psi became non-pos.def. during iterations")
	else if(msg==6)
	  warning("Value of sigma became non-pos.def. during iterations")
        else if(msg==7)
          warning("log-density not concave at one or more scoring steps")}
        llvec<-tmp$llvec[1:iter]
	llovec<-tmp$llovec[1:iter]
	converged<-cvgd==as.integer(1)
	if(!converged) warning(paste("did not converge by",
	   format(iter),"iterations"))
	#
	list(beta=beta,sigma=sigma,psi=psi,eb=tmp$eb,varb=tmp$varb,xtwxinv=tmp$xtwxinv,
			converged=converged,iter=iter,npatt=npatt,pstfin=pstfin,iposn=iposn,patt=patt,rmat=rmat,
			logll=llvec,logoll=llovec,clock=clock)}
