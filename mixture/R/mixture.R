.packageName<-'mixture'


npar.model <- function(modelname=NULL, p=NULL, G=NULL) {
	val = numeric(3)
	val[1] = G-1 
	val[2] = G*p
	val[3] = ncovpar(modelname= modelname, p=p, G=G)
	val = sum(val)
	return(val)
}

ncovpar <- function(modelname=NULL, p=NULL, G=NULL) {
	if (is.null(p)) stop("p is null")
	if (is.null(G)) stop("G is null")
	if (is.null(modelname)) stop("modelname is null")

	     if (modelname == "EII") npar = 1
	else if (modelname == "VII") npar = G
	else if (modelname == "EEI") npar = p
	else if (modelname == "VEI") npar = p + G -1	
	else if (modelname == "EVI") npar = p*G - G +1
	else if (modelname == "VVI") npar = p*G
	else if (modelname == "EEE") npar = p*(p+1)/2
	else if (modelname == "EEV") npar = G*p*(p+1)/2 - (G-1)*p	
	else if (modelname == "VEV") npar = G*p*(p+1)/2 - (G-1)*(p-1)
	else if (modelname == "VVV") npar = G*p*(p+1)/2
	else if (modelname == "EVE") npar = p*(p+1)/2 + (G-1)*(p-1)
	else if (modelname == "VVE") npar = p*(p+1)/2 + (G-1)*p
	else if (modelname == "VEE") npar = p*(p+1)/2 + (G-1)
	else if (modelname == "EVV") npar = G*p*(p+1)/2 - (G-1)
	else stop("modelname is not correctly defined")
	
	return(npar)		
}

model.type <- function(modelname=NULL, Sk=NULL, ng=NULL, D=NULL, mtol=1e-10, mmax=10) {
	if (is.null(modelname)) stop("modelname is null")

	     if (modelname == "EII") val = msEII(Sk=Sk, ng=ng)
	else if (modelname == "VII") val = msVII(Sk=Sk, ng=ng)
	else if (modelname == "EEI") val = msEEI(Sk=Sk, ng=ng)
	else if (modelname == "VEI") val = msVEI(Sk=Sk, ng=ng, eplison= mtol, max.iter= mmax)
	else if (modelname == "EVI") val = msEVI(Sk=Sk, ng=ng)
	else if (modelname == "VVI") val = msVVI(Sk=Sk, ng=ng)
	else if (modelname == "EEE") val = msEEE(Sk=Sk, ng=ng)
	else if (modelname == "EEV") val = msEEV(Sk=Sk, ng=ng)
	else if (modelname == "VEV") val = msVEV(Sk=Sk, ng=ng, eplison= mtol, max.iter= mmax)
	else if (modelname == "VVV") val = msVVV(Sk=Sk, ng=ng)
	else if (modelname == "EVE") val = msEVE(Sk=Sk, ng=ng, D0=D, eplison= mtol, max.iter= mmax)
	else if (modelname == "VVE") val = msVVE(Sk=Sk, ng=ng, D0=D, eplison= mtol, max.iter= mmax)
	else if (modelname == "VEE") val = msVEE(Sk=Sk, ng=ng, eplison= mtol, max.iter= mmax)
	else if (modelname == "EVV") val = msEVV(Sk=Sk, ng=ng)
	else stop("modelname or covtype is not correctly defined")
	
	if (!is.list(val)) val = list(sigma=val)
	return(val)		
}


gpcm <- function(data=NULL,  G=1:3, mnames=NULL, start=0, label=NULL, veo=FALSE, nmax=1000, atol=1e-8, mtol=1e-8, mmax=10, pprogress=FALSE, pwarning=FALSE) {
        set.seed(102)
	if (is.null(data)) stop('Hey, we need some data, please! data is null')
	if (!is.matrix(data)) stop('The data needs to be in matrix form')
	if (!is.numeric(data)) stop('The data is required to be numeric')
	if (nrow(data) == 1) stop('nrow(data) is equal to 1')
	if (ncol(data) == 1) stop('ncol(data) is equal to 1; This function currently only works with multivariate data p > 1')
	if (any(is.na(data))) stop('No NAs allowed.')
	if (is.null(G)) stop('G is NULL')
	G = as.integer(ceiling(G))
	if (!is.integer(G)) stop('G is not a integer')
	if ( any(G < 1)) stop('G is not a positive integer')
	
	if (is.null(mnames) )  mnames = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", "VEV", "VVV", "EVE", "VVE", "VEE", "EVV")
#	if (is.null(mnames) )  mnames = c("EVE", "VVE", "VEE", "EVV")

	bic = array(0, dim= c(length(G), length(mnames), 3), dimnames=list(G, mnames, c('loglik', "npar", "BIC")) )
	#BIC = matrix(0, nrow=length(G), ncol=length(mnmaes), dimnames=list(G, mnames) )	
	model = NULL; curBIC = Inf;
	for (g in 1:length(G)) {
	for (i in 1:length(mnames)) {
		if ( pprogress ) print(c(G[g],mnames[i]))
		if (veo | npar.model(modelname=mnames[i], p=ncol(data), G=G[g]) < nrow(data)) {
        
 #  print("gpcm EM1")
        a =	try( { EM(data=data, G=G[g], nmax=nmax, covtype=mnames[i], start=start, label=label, atol= atol, mtol=mtol, mmax=mmax ) }, silent= !pwarning)
 #  print("gpcm EM2")
#        a =	EM(data=data, G=G[g], nmax=nmax, covtype=mnames[i], start=start, label=label, atol= atol, mtol=mtol, mmax=mmax ) 
		if(length(a) > 1){
			bic[g,i,1:2] = c(a$loglikn, npar.model(modelname=mnames[i], p=ncol(data), G=G[g]) )
			bic[g,i,3] =  -2*bic[g,i,1] + bic[g,i,2]*log(nrow(data)) 
#			cat(bic[g,i,1],bic[g,i,3], curBIC, "\n")
			if (is.nan(bic[g,i,3]) | is.infinite(bic[g,i,3]) ){
				bic[g,i,1:2] = c(NA, npar.model(modelname=mnames[i], p=ncol(data), G=G[g]) )
				bic[g,i,3] =  NA
			}
			else if ( bic[g,i,3] < curBIC) {
					model  = append(a,list(mtype=mnames[i]))
					curBIC = bic[g,i,3] 
			}
		} else {
			bic[g,i,1:2] = c(NA, npar.model(modelname=mnames[i], p=ncol(data), G=G[g]) )
			bic[g,i,3] =  NA
		}
	} else {
		bic[g,i,1:2] = c(NA, npar.model(modelname=mnames[i], p=ncol(data), G=G[g]) )
		bic[g,i,3] =  NA	
	} 	
	
	}}
	
	if ( is.null(start) ) startobject= "deterministic annealing with default values."
	else if ( is.matrix(start)) startobject= "a user specified initialization  matrix."
	else if ( is.function(start)) startobject= "a user specified initialization function."
	else if ( length(start) > 1 ) startobject= "deterministic annealing with user specified values"
	else if ( start == 0) startobject= "k-means"
	else if ( start > 0) startobject= paste( start, " random initializations", collapse="" )
	
	bicModel = list(G=model$G, covtype=model$mtype, bic=curBIC )
	val = list( start=start, startobject= startobject, gpar=model$gpar, loglik=model$loglik, z=model$z, map=model$map, BIC=bic, bicModel= bicModel)
#	val = list( model=model, BIC=bic, bicModel= bicModel)
	
    class(val)<-"gpcm"
	return(val)
}


print.gpcm <-function(x, ...){
#    a = object$BIC[,,1]
    cat("The model choosen by applying the BIC criteria has", x$bicModel$G, "component(s) and the", x$bicModel$covtype, "covariance structure\n using", x$startobject, "starts. \n"  )
    #endPrint(a)
}
summary.gpcm <- function(object, ...){
    bicl = object$BIC[,,3]
    cat("BIC for each model, number of components (rows), and covariance structure (columns).\n")
    print.default( bicl)
}

plot.gpcm <- function(x, ...) {
	bicl = x$BIC[,,1]
	g   = dimnames(bicl)[[1]]
	cov = dimnames(bicl)[[2]]
	ncov = length(cov)
	ng   = length(g)

	plot( g, bicl[,1],  ylim=range(bicl,na.rm=TRUE), xlab="# components", ylab="BIC",type='l')
	if (ncov >1) {
		for (i in 2:ncov) lines( g, bicl[,i], col=i)
	}
	legend( "bottomright", legend=cov,  col=1:ncov,lty=1)
	
}



igpar_1 <- function(data=NULL, g=NULL, covtype=NULL, start=NULL, labels=NULL, mtol=NULL, mmax=NULL) {
        if (is.null(start)) start = seq( 0.1, 1, length.out=25)

        if (is.function(start) ) w = start(data=data, g=g, covtype=covtype)
        else if ( !is.null(dim(start)) ) w = start
        else if (length(start) > 1) w = igparv(data=data, g=g, covtype=covtype, vseq0=start , mtol=mtol, mmax = mmax, labels=labels)
        else if (start == 0) w = igpark(data=data, g=g, covtype=covtype)
        else if (start >  0) w = rgpar(data=data, g=g, covtype=covtype, n=ceiling(start), labels=labels, mtol=mtol, mmax = mmax )
        else stop(paste('Initialization method ', start, " does not compute!"))

        if (!is.matrix(w)) stop("The zij initialization matrix is not a matrix")
        if (any(w < 0)) stop("Some of the elements of the zij initialization matrix are less than zero")
#       if ( any(apply(w,1,sum) !=1)  ) stop("Some of the rows of the zij initialization matrix do not sume to 1")
        if (nrow(w) < nrow(data) ) stop("The nrow(zij) of the initialization matrix is less than nrow(data)")
        if (nrow(w) > nrow(data) ) stop("The nrow(zij) of the initialization matrix is greater than nrow(data)")
        if (ncol(w) < g ) stop("The nrow(zij) of the initialization matrix is less than g")
        if (ncol(w) > g ) stop("The nrow(zij) of the initialization matrix is greater than g")
        w = combinewk(w, label=labels)

        return(w)
}



igparv <- function(data=NULL, g=NULL, covtype=NULL, vseq0=NULL, labels=NULL, mtol=NULL, mmax = NULL) {
	vseq0 = as.numeric(vseq0)
	if (is.null(vseq0)) stop('The sequence for deterministic annealing is NULL')
	if ( !all( vseq0 <=1 &  vseq0 >=0)  ) stop('The sequence for deterministic annealing must be between 0 and 1')
	w = EMv(data=data, G=g, vseq=vseq0, m=1, label=labels, covtype = covtype, mtol= mtol, mmax = mmax )
	return(w)
}



rgpar <- function(data=NULL, g=NULL, covtype=NULL, n=1, labels=NULL, mtol=NULL, mmax = NULL) {
	rn = 5;
	w.old    = rwgpar(data=data, g=g,covtype= covtype, labels=labels)
	temp.old = EMn(data=data, G=g, w=w.old, label= labels, covtype= covtype, n=rn, atol=1e-14, mtol= mtol, mmax= mmax)

	for (i in 1:n) { 
		w.new    = rwgpar(data=data, g=g,covtype= covtype, labels=labels)
		temp.new = EMn(data=data, G=g, w=w.new, label= labels, covtype= covtype, n=rn, atol=1e-14, mtol= mtol, mmax= mmax)
		if (temp.new$loglikn > temp.old$loglikn) temp.old = temp.new
	}	

	return(temp.old$z)
}



rwgpar <- function(data=NULL, g=NULL,covtype=NULL, labels=NULL) {
	if (g > 1) {
		w = matrix(rexp(nrow(data)*g),nrow=nrow(data),ncol=g)
		sw = apply(w, 1, sum)
		w  = sweep(w, 1, sw, '/')
		#w = matrix(t(apply(w,1, function(z) { z/sum(z)})),nrow=nrow(data),ncol=g)
		w = combinewk(w, label=labels)
	} else {
		w = matrix(1,nrow=nrow(data),ncol=g) 
	}
	return(w)
}


igpark <- function(data=NULL, g=NULL, covtype=NULL) {
	if (g==1) w = matrix(0,nrow=nrow(data),ncol=1)
	else {
	lw = kmeans(x=data, centers=g, iter.max= 25 )$cluster
	w  = combinewk(matrix(0,nrow=nrow(data),ncol=g), label=lw)
	}
	return(w)
}





run_em <- function(x=NULL, G=NULL, z=NULL, nmax=NULL, atol=NULL, mtol=NULL, mmax=NULL, label=NULL, covtype=NULL){
	if (!is.matrix(data)) as.matrix(x)
    N = nrow(x)
    p = ncol(x)
    if (is.null(label) ) label = rep(1,N)	
    MAPP = numeric(N)
    logl = numeric(nmax)
    counter = 0
    sigmar    = matrix(0, nrow=G, ncol=p^2 )
    invsigmar = matrix(0, nrow=G, ncol=p^2 )
    mu = matrix(0, nrow = G, ncol=p )    
    D = matrix(0, nrow = p, ncol=p )    
    pi = numeric(G)
    
    temp_em<-.C("main_loop", as.integer(N), as.integer(p), as.integer(G), as.double(z),
        as.double(sigmar), as.double(invsigmar), as.double(mu), as.double(pi), 
        as.integer(nmax), as.double(atol), as.double(mtol),
        as.integer(mmax), as.double(x), as.integer(label), as.character(covtype), 
        as.double(logl), as.integer(counter), as.integer(MAPP), as.double(D), PACKAGE="mixture")

    z        = matrix(temp_em[[4]], nrow=N, ncol=G)
	num.iter = temp_em[[17]]
    map      = temp_em[[18]]
    DD       = temp_em[[19]]
    loglik   = temp_em[[16]]
    loglik   = loglik[1:num.iter]
    loglikn  = loglik[num.iter]

    sigma    = array(temp_em[[5]], dim= c(p,p,G) ) 
    invsigma = array(temp_em[[6]], dim= c(p,p,G) )
    mu       = matrix(temp_em[[7]], nrow=p, ncol=G, byrow=TRUE)
    pi       = temp_em[[8]]

    gpar= list()
	for (k in 1:G ) {
		gpar[[k]] = list()		
		gpar[[k]]$mu       = mu[,k]
		gpar[[k]]$sigma    = sigma[,,k]
		gpar[[k]]$invSigma = invsigma[,,k]
		gpar[[k]]$logdet   = log(det(sigma[,,k]))

	}
	gpar$pi = temp_em[[8]]	
	if (covtype == "EVE") gpar$D  = matrix(DD, nrow=p, ncol=p)
	if (covtype == "VVE") gpar$D  = matrix(DD, nrow=p, ncol=p)	

    val = list(z=z, loglik=loglik, gpar = gpar, loglikn=loglikn, num.iter = num.iter, map=map, G=G, mtype= covtype)
    return(val)
}




EMp1 <- function(data=NULL, covtype=NULL) {

	d= ncol(data);
	w = rep(1, nrow(data))
	temp = cov.wt(data, wt=w, center=TRUE, method="ML")
	
	gpar= list()	
	gpar[[1]]          = list()
	gpar[[1]]$mu       = temp$center
	if (substr(covtype,2,2) == 'I') {
		c0 = mean(diag(temp$cov))
		gpar[[1]]$sigma    = c0*diag(d)
		gpar[[1]]$invSigma = 1/c0*diag(d)
		gpar[[1]]$logdet   = d*log( c0 ) 
		
	} else if (substr(covtype,3,3) == 'I') {
		c0 = diag(temp$cov)
		gpar[[1]]$sigma    = diag(c0)   
		gpar[[1]]$invSigma = diag(1/c0)
		gpar[[1]]$logdet   = sum( log( c0 ) )

	} else {
		gpar[[1]]$sigma    = temp$cov
		gpar[[1]]$invSigma = solve(gpar[[1]]$sigma)
		gpar[[1]]$logdet   = log(det(gpar[[1]]$sigma))
	}	
	gpar$pi = 1

	zlog = -1/2*mahalanobis(x=data, center=gpar[[1]]$mu, cov=gpar[[1]]$invSigma, inverted=TRUE) -1/2*gpar[[1]]$logdet - d/2*(log(2)+log(pi))
	loglik = sum(zlog)	
	
	val = list(gpar=gpar, z=matrix(w,nrow=nrow(data), ncol=1), iterations=1, loglik=loglik, loglikn=loglik, map=w, G=1, mtype=covtype)
	return(val)
	}



EM <- function(data=NULL, G=2, start=1, label=NULL, covtype=NULL, nmax=1000, atol=1e-8, mtol=1e-8, mmax=10 ) {
	if (G == 1) {
		val = EMp1(data=data, covtype=covtype) 	
	} else {	
		# G>1
    	w   = igpar_1(data=data, g=G, covtype=covtype, start=start, labels=label, mtol=mtol, mmax=mmax)
	    val = run_em(x=data, G=G, z=w, nmax=nmax, atol=atol, mtol=mtol, mmax=mmax, label=label,covtype=covtype)
     } 
     
     return(val)
}



EMn <- function(data=NULL, G=2, w=NULL, label=NULL, covtype=NULL, n=5, atol=1e-14, mtol=1e-8, mmax=10 ) {
	 nmax = n 
     if (is.null(w)) stop("EMn has a w that is null") 
     val <- run_em(G=G,z=w,nmax=nmax,atol=atol, mtol =mtol,mmax=mmax,x=data,label=label,covtype=covtype)
     return(val)
}


EMv <- function(data=NULL, G=3, vseq=c(1,1), m=2, label=NULL, covtype="VVV", mtol=NULL, mmax = NULL  ) {
 	tempw = rwgpar(data=data, g=G, covtype= covtype, labels=label)
 	gpar  = m.step(data=data, covtype=covtype, w=tempw, D=NULL, mtol= mtol, mmax = mmax)
 	tempw = e.step(data=data, gpar=gpar, labels=label, v=vseq[1])
	llik  = numeric(length(vseq)*m)
	for (i in 1:length(vseq)) { for (j in 1:m) {		
 		gpar    = m.step(data=data, covtype=covtype, w=tempw, D=gpar$D, mtol= mtol, mmax = mmax)
 		tempw   = e.step(data=data, gpar=gpar, labels=label, v=vseq[i])
	}}
	return(tempw)
	}



weights <- function(data=NULL, gpar=NULL, v=1) {
	d = ncol(data)
	G = length(gpar$pi)	
	if (G > 1) {
		zlog = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
		for (k in 1:G ) zlog[,k] = -1/2*mahalanobis(x=data, center=gpar[[k]]$mu, cov=gpar[[k]]$invSigma, inverted=TRUE) -1/2*gpar[[k]]$logdet - d/2 *log(2*pi)

		w = t(apply( zlog, 1, function(z,wt,v) { 
			x= exp( v*(z + log(wt)) ) 
			x=x/sum(x);
			return(x) }, wt=gpar$pi,v=v ))
	} else w = matrix(1,nrow=nrow(data), ncol=G)
	return(w)
	}


e.step <- function(data=NULL, gpar=NULL, labels=NULL, v=1) {
		w = weights(data=data, gpar=gpar,v=v)
		if (!is.null(labels)) w = combinewk(weights=w, label= labels)
		return(w)
}


m.step <- function(data=NULL, covtype=NULL, w=NULL, D=NULL, mtol=NULL, mmax=NULL) {
	G= ncol(w);
	d= ncol(data);
	Sk = array(0, c(d,d,G) )
	gpar= list()	
	for (k in 1:G ) {
		gpar[[k]] = list()		
		temp = cov.wt(data, wt=w[,k], center=TRUE, method="ML")
		gpar[[k]]$mu    = temp$center
		if (!any(is.na( temp$cov))) gpar[[k]]$sigma = temp$cov
		Sk[,,k]   = temp$cov
	}
	gpar$pi = apply(w,2,mean)

	temp = model.type(modelname = covtype, Sk=Sk, ng=gpar$pi, D=D, mtol= mtol, mmax= mmax )
	gpar$D = temp$D
	for (k in 1:G ) {
		gpar[[k]]$sigma    = temp$sigma[,,k] #+ diag( 1-v , d, d)
		gpar[[k]]$invSigma = temp$invSigma[,,k] 
		gpar[[k]]$logdet   = temp$logdet[k]
		}
	return(gpar)
	}


	
loglik <- function(data, gpar) {
	# output is a G x nrow(data) matrix
	d = ncol(data)
	G = length(gpar$pi)
	zlog = matrix(0, nrow=nrow(data), ncol=G)
	for (k in 1:G) zlog[,k] = -1/2*mahalanobis(x=data, center=gpar[[k]]$mu, cov=gpar[[k]]$invSigma, inverted=TRUE) -1/2*gpar[[k]]$logdet - d/2*(log(2)+log(pi))
	
	w = apply( exp(zlog),1,function(z,wt) { sum(z*wt) } , wt=gpar$pi)
	val = sum(log(w))
	if( is.nan(val) ) val = NA	
	return(val)
	}
	

MAP <- function(data, gpar, label=NULL) {
	w = weights(data=data, gpar=gpar, v=1)
	if (!is.null(label)) w = combinewk(weights=w, label= label)
	z = apply(w, 1, function(z) { z=(1:length(z))[z==max(z)]; return(z[1]) })
	z = as.numeric(z)
	return( z)	
	}
	

combinewk <- function(weights=NULL, label=NULL)	{
#        cat("label is",label,"\n")
	# known is a numeric with 
	# 0 if unknown group membership 
	# 1,2,3,.. for label of known group
	if (!is.null(label)) { #stop('label is null')
	label = as.integer(label)
	
	if (any(!is.integer(label)))	stop("Labels are not integers")
	if (any(label < 0 ))	 stop("Labels can only be positive integers")
	if (ncol(weights) < max(label) ) stop("Number of groups is less then the number groups given by labels")
			
	if ( sum(label!=0) == nrow(weights) ) {
		if (ncol(weights) > max(label) ) stop("Every observations has a label; Cannot fit more groups to the data then given the by the labels.")
	}

	kw     = label !=0
 	for (j in 1:ncol(weights)) weights[kw,j] = (label == j)[kw]
 	}
	return(weights)	
}






msEEE <- function(Sk=NULL, ng=NULL) {
	# Sk is an array of with dim (p x p x G)
	# ng is a vector of length G of weights of Wk
	d = dim(Sk)[1];
	G = dim(Sk)[3];
	W = sumSk.wt(Sk=Sk, wt=ng, d=d, G=G)/sum(ng)

	val = array(0, c(d,d,G))
	for (g in 1:G) val[,,g] = W

	logdetW = log(det(W))
	invW    = solve(W)
	val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
	for (g in 1:G) { 
		val$sigma[,,g]    = W
		val$invSigma[,,g] = invW
		val$logdet[g]     = logdetW
		}
	return(val)
}

sumSk.wt <- function(Sk=NULL, wt=NULL, d=NULL, G=NULL) {
	# Sum Sk over the groups with weights wt.
	W = matrix(0, nrow=d, ncol=d )
	for (g in 1:G) W = W + Sk[,,g]* wt[g]
	return(W)
}


msEEV <- function(Sk=NULL, ng=NULL, eplison=1e-20, max.iter= 100) {
	# Wk is a list of length G with matrix (p x p)
	# ng is a vector of length G of weights of Wk
	d = dim(Sk)[1];
	G = dim(Sk)[3];
	
	EWk = Sk
	A  = matrix(0,d,d)
	for (g in 1:G) {
		Wk = Sk[,,g]*ng[g]
		EWk[,,g] = eigen(Wk)$vectors
		A = A + t(EWk[,,g]) %*% Wk %*% EWk[,,g]
	}
	
	lam = prod(diag(A))^(1/d)
	A   = A/lam
	lam = lam/sum(ng)

	val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
	for (g in 1:G) { 
		val$sigma[,,g]    =   lam *( EWk[,,g] %*% A %*% t(EWk[,,g]) )
		val$invSigma[,,g] = 1/lam *( EWk[,,g] %*% diag(1/diag(A),d) %*% t(EWk[,,g]) )
		val$logdet[g]     =  d*log(lam) 
		}
	return(val)
}



getA <- function(Ok=NULL, lam=NULL, G=NULL, d=NULL) {
	A  = matrix(0,d,d)
	for (g in 1:G) A = A + Ok[,,g]/lam[g]
	A= diag(A)
	A= diag( A/prod(A)^(1/d) )
	return( A )	
}


getOk <- function(Sk=NULL, ng=NULL, G=NULL) {
	Ok  = Sk
	for (g in 1:G) {
		Wk  = Sk[,,g]*ng[g]
		EWk = eigen(Wk)$vectors
		Ok[,,g] = t(EWk) %*% Wk %*% EWk
		Ok[,,g] = diag(diag(Ok[,,g]))
	}
	return(Ok)	
}

getEkOk <- function(Sk=NULL, ng=NULL, G=NULL) {
	Ok = Sk
	EWk = Sk
	for (g in 1:G) {
		Wk  = Sk[,,g]*ng[g]
		EWk[,,g] = eigen(Wk)$vectors
		Ok[,,g] = t(EWk[,,g]) %*% Wk %*% EWk[,,g]
	}
	return(list(Ok=Ok,EWk=EWk))	
}


msVEV <- function(Sk=NULL, ng=NULL, eplison=1e-14, max.iter= 100) {
	# Wk is a list of length G with matrix (p x p)
	# ng is a vector of length G of weights of Wk
	d = dim(Sk)[1];
	G = dim(Sk)[3];
	
	temp = getEkOk(Sk=Sk, ng=ng, G=G)
	Ok  = temp$Ok
	EWk = temp$EWk
	lam = apply(Ok,3,function(z) { sum(diag(z)) } )/(ng*d)
	A   = getA(Ok=Ok, lam=lam, d=d, G=G) 
	lam = apply(Ok,3,function(z, invA) { sum(diag(z*invA)) }, invA=diag(1/diag(A)) )/(ng*d)

	conv = c( d*sum(ng*(1+log(lam))), Inf  )
	count = 1 
	while ( diff(conv)/conv[1] > eplison & count < max.iter) {
		A   = getA(Ok=Ok, lam=lam, d=d, G=G) 
		lam = apply(Ok,3,function(z, invA) { sum(diag(z*invA)) }, invA=diag(1/diag(A)) )/(ng*d)
		conv = c(d*sum(ng*(1+log(lam))), conv[1] )
		count = count +1
	}

	val = array(0, c(d,d,G))
	for (g in 1:G) val[,,g] = lam[g] * ( EWk[,,g] %*% A %*% t(EWk[,,g]) )
	
	val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
	for (g in 1:G) { 
		val$sigma[,,g]    =    lam[g] * ( EWk[,,g] %*% A %*% t(EWk[,,g]) )
		val$invSigma[,,g] =  1/lam[g] * ( EWk[,,g] %*% diag(1/diag(A),d) %*% t(EWk[,,g]) )
		val$logdet[g]     =  d*log(lam[g])
		}
	return(val)
}



msVVV <- function(Sk=NULL, ng=NULL) {
	# Wk is a list of length G with matrix (p x p)
	# ng is a vector of length G of weights of Wk
	d = dim(Sk)[1];
	G = dim(Sk)[3];

	val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
	for (g in 1:G) { 
		val$sigma[,,g]    =  Sk[,,g]
		val$invSigma[,,g] =  solve(Sk[,,g])
 		val$logdet[g]     =  log(det(Sk[,,g]) )
		}
	return(val)
}


msEEI <- function(Sk=NULL, ng=NULL) {
	# Wk is a list of length G with matrix (p x p)
	# ng is a vector of length G of weights of Wk
	d = dim(Sk)[1];
	G = dim(Sk)[3];
	
	W = sumSk.wt(Sk=Sk, wt=ng, d=d, G=G)/sum(ng)
	B = diag(diag(W))

	val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G) )
	for (g in 1:G) { 
		val$sigma[,,g]    = B
		val$invSigma[,,g] = diag( 1/diag(B),d )
		val$logdet[g] =  sum(log( diag(B) ))
		}
	return(val)
}





msVEI <- function(Sk=NULL, ng=NULL, eplison=1e-20, max.iter= 100) {
	# Wk is a list of length G with matrix (p x p)
	# ng is a vector of length G of weights of Wk
	d = dim(Sk)[1];
	G = dim(Sk)[3];
	
	lam = apply(Sk,3,function(z) { sum(diag(z)) } )/d
#	W = sumSk.wt(Sk=Sk, wt=ng*lam, d=d, G=G)
	W = sumSk.wt(Sk=Sk, wt=ng/lam, d=d, G=G)
	W = diag(W)
	B = diag(W/prod(W)^(1/d))
	lam = apply(Sk,3,function(z, invB) { sum(diag(z*invB)) }, invB=diag(1/diag(B) ) )/d

	conv = c( d*sum(ng*(1+log(lam))), Inf  )
	count = 1 
	while ( abs(diff(conv)) > eplison & count < max.iter) {
#	while ( count < max.iter) {
#		W = sumSk.wt(Sk=Sk, wt=ng*lam, d=d, G=G)
		W = sumSk.wt(Sk=Sk, wt=ng/lam, d=d, G=G)

		W = diag(W)
		B = diag(W/prod(W)^(1/d))
		lam = apply(Sk,3,function(z, invB) { sum(diag(z*invB)) }, invB=diag(1/diag(B) ) )/d

		conv = c(d*sum(ng*(1+log(lam))), conv[1] )
		count = count +1
	}

	val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
	for (g in 1:G) { 
		val$sigma[,,g]    = lam[g]*B
		val$invSigma[,,g] = diag( 1/diag(B) * 1/lam[g],d )
		val$logdet[g]     = d*log(lam[g])
		}
	return(val)
}


msEVI <- function(Sk=NULL, ng=NULL) {
	# Wk is a list of length G with matrix (p x p)
	# ng is a vector of length G of weights of Wk
	d = dim(Sk)[1];
	G = dim(Sk)[3];
	
	Bk = matrix(0,nrow=d, ncol=G)
	for (g in 1:G) Bk[,g] = diag(Sk[,,g]*ng[g])
	lam = apply(Bk, 2, prod)^(1/d) 
	Bk  = sweep(Bk, 2, 1/lam, FUN="*")	
	lam = sum(lam)/sum(ng)

	val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
	for (g in 1:G) { 
		val$sigma[,,g]    = lam * diag(Bk[,g],d)
		val$invSigma[,,g] = 1/lam * diag(1/Bk[,g],d)
		val$logdet[g]     = d*log(lam)
		}
	return(val)
}




msVVI <- function(Sk=NULL, ng=NULL) {
	# Wk is a list of length G with matrix (p x p)
	# ng is a vector of length G of weights of Wk
	d = dim(Sk)[1];
	G = dim(Sk)[3];
	
#	Bk = matrix(0,nrow=d, ncol=G)
#	for (g in 1:G) Bk[,g] = diag(Sk[,,g]*ng[g])
	
	val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
	for (g in 1:G) { 
		val$sigma[,,g]    = diag(diag(Sk[,,g]),d)
		val$invSigma[,,g] = diag(1/diag(Sk[,,g]),d)
		val$logdet[g]     = sum( log(diag(Sk[,,g])) )
		}
	return(val)
}

msEII <- function(Sk=NULL, ng=NULL) {
	# Wk is a list of length G with matrix (p x p)
	# ng is a vector of length G of weights of Wk
	d = dim(Sk)[1];
	G = dim(Sk)[3];
	
	W = sumSk.wt(Sk=Sk, wt=ng, d=d, G=G)/sum(ng)
	lam = sum(diag(W))/(sum(ng)*d)
	
	val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G) )
	for (g in 1:G) { 
		val$sigma[,,g] = diag(rep(lam,d), d)
		val$invSigma[,,g] = diag(rep(1/lam,d),d)
		val$logdet[g] = d*log(lam) 
		}
	return(val)
}


msVII <- function(Sk=NULL, ng=NULL) {
	# Wk is a list of length G with matrix (p x p)
	# ng is a vector of length G of weights of Wk
	d = dim(Sk)[1];
	G = dim(Sk)[3];
	
	val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
	for (g in 1:G) { 
		sumdiagSkg = sum(diag(Sk[,,g]))
		val$sigma[,,g] = diag(rep(sumdiagSkg/d,d))
		val$invSigma[,,g] = diag(rep( d/sumdiagSkg, d))
		val$logdet[g] =  d* log( sumdiagSkg ) - d*log(d)
		}
	return(val)
}








msVVE <- function(Sk=NULL, ng=NULL) {
	# Sk is a list of length G with matrix (p x p)
	# ng is a vector of length G of weights of Wk
	d = dim(Sk)[1];
	G = dim(Sk)[3];
	
#	lam = numeric(G)
#	for (g in 1:G) 	lam[g] = sum(diag(Sk[,,g]))/d
	
	val = array(0, c(d,d,G))
	for (g in 1:G) val[,,g] = diag(rep(sum(diag(Sk[,,g]))/d,d))
	return(val)
}



##################### 
### new.cov models
##################### 



newD3.MM <- function(D=NULL, d=NULL, G=NULL, Wk=NULL, Ak=NULL, tmax = 100) {
	z = matrix(0,d,d)
	lambda =0 
	for (g in 1:G) {
		lambdak = max(eigen(Wk[,,g])$values)
		z = z + diag(1/Ak[,g]) %*% t(D) %*% Wk[,,g]  -lambdak *(diag(1/Ak[,g])%*% t(D) )
	} 
	z1 = svd(z)
	Xk1 = (z1$v) %*% t(z1$u) 
	return( Xk1 )
}

newD4.MM <- function(D=NULL, d=NULL, G=NULL, Wk=NULL, Ak=NULL, tmax = 100) {
	z = matrix(0,d,d)
	lambda =0 
	for (g in 1:G) {
		lambdak = max(1/Ak[,g])
		#z = z + diag(1/Ak[,g]) %*% t(D) %*% Wk[,,g]  -lambdak *(diag(1/Ak[,g])%*% t(D) )
		z = z + Wk[,,g] %*% (D) %*% diag(1/Ak[,g])  -lambdak *(Wk[,,g]%*% (D) )
	} 
	z1 = svd(z)
#	Xk1 = (z1$v) %*% t(z1$u) 
#	return( t(Xk1) )
# OR 
	Xk1 = (z1$v) %*% t(z1$u) 
	return( t(Xk1) )

}

newD <- function(D=NULL, d=NULL, G=NULL, Wk=NULL, Ak=NULL, tmax = 100) {
	D6 =D
	D6 = newD3.MM(D=D6, d=d, G=G, Wk=Wk, Ak=Ak, tmax = 100)
	D6 = newD4.MM(D=D6, d=d, G=G, Wk=Wk, Ak=Ak, tmax = 100)
	return(D6)
}



testval <- function(Wk=NULL, Ak=NULL, D=NULL, G=NULL) {
	z = numeric(G)
#	for (g in 1:G) z[g] = sum(diag( D %*%  diag(1/Ak[,g]) %*% t(D) %*% Wk[,,g]))
	for (g in 1:G) z[g] = sum(diag( t(D) %*% Wk[,,g] %*% D %*%  diag(1/Ak[,g])  ))
	return(sum(z))
}



testgrad.D <- function(D=NULL, d=NULL, G=NULL, Wk=NULL, Ak=NULL) {
	z = matrix(0,d,d)
	for (g in 1:G) z = z + Wk[,,g] %*% D %*% diag(1/Ak[,g]) 
	return(2*z)
}



msEVE <- function(Sk=NULL, ng=NULL, eplison=1e-20, max.iter= 10, D0=NULL) {
	# Wk is a list of length G with matrix (p x p)
	# ng is a vector of length G of weights of Wk
	d = dim(Sk)[1]; G = dim(Sk)[3];

	Wk = Sk
	W  = matrix(0,d,d)
	for (g in 1:G) {
		Wk[,,g] = Sk[,,g]*ng[g]
		W  = W + Wk[,,g]
	}
#	D  = diag(rep(1,d))
	D = D0
	if (is.null(D)) D  = diag(rep(1,d))
#	D  = t(eigen(W)$vectors)
	Ak = apply(Wk,3, function(z,D) { diag( t(D) %*% z %*% (D) ) }, D=D )
	Ak = apply(Ak,2,function(z) { z/prod(z)^(1/length(z))})
#print(c(0, 1, testval(Wk=Wk,Ak=Ak,D=D,G=G)))		
	D = newD(D=D, Wk=Wk, Ak=Ak, G=G, d=d)
#print(c(0, 2, testval(Wk=Wk,Ak=Ak,D=D,G=G)))		

 	conv = c( testval(Wk=Wk,Ak=Ak,D=D,G=G), Inf  )
	count = 1  
	while ( diff(conv)/abs(conv[1]) > eplison & count < max.iter) {
		D = newD(D=D, Wk=Wk, Ak=Ak, G=G, d=d) 
		Ak = apply(Wk,3, function(z,D) { diag( t(D) %*% z %*% (D) ) }, D=D )
		Ak = apply(Ak,2,function(z) { z/prod(z)^(1/length(z))})
	
#print(c(count, 0, testval(Wk=Wk,Ak=Ak,D=D,G=G)))		
		conv = c(testval(Wk=Wk,Ak=Ak,D=D,G=G), conv[1] )
		count = count +1
	}
	lam =  0
	for (g in 1:G) lam = lam  + sum(diag( D %*% diag(1/Ak[,g])%*% t(D) %*% Wk[,,g] ))
	lam = lam/(sum(ng)*d)

	val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G), D=D  )
	for (g in 1:G) { 
		val$sigma[,,g]    = D %*% diag(lam*Ak[,g])%*% t(D)
		val$invSigma[,,g] = D %*% diag(1/lam*1/Ak[,g])%*% t(D)
		val$logdet[g]     = d*log( lam ) 
		}
	return(val)

}


msVVE <- function(Sk=NULL, ng=NULL, eplison=1e-20, max.iter= 10, D0=NULL) {
	# Wk is a list of length G with matrix (p x p)
	# ng is a vector of length G of weights of Wk
	d = dim(Sk)[1]; G = dim(Sk)[3];

	Wk = Sk
	W  = matrix(0,d,d)
	for (g in 1:G) {
		Wk[,,g] = Sk[,,g]*ng[g]
		W  = W + Wk[,,g]
	}

	D = D0
	if (is.null(D)) D  = diag(rep(1,d))

	Ak = apply(Wk,3, function(z,D) { diag( t(D) %*% z %*% D ) }, D=D )
	Ak = apply(Ak,2,function(z) { z/prod(z)^(1/length(z))})

#print(c(0, 1, testval(Wk=Wk,Ak=Ak,D=D,G=G)))		
	D = newD(D=D, Wk=Wk, Ak=Ak, G=G, d=d)
#print(c(0, 2, testval(Wk=Wk,Ak=Ak,D=D,G=G)))		

 	conv = c( testval(Wk=Wk,Ak=Ak,D=D,G=G), Inf  )
	count = 1  
	while ( diff(conv)/abs(conv[1]) > eplison & count < max.iter) {
		Ak = apply(Wk,3, function(z,D) { diag( t(D) %*% z %*% D ) }, D=D )
		Ak = apply(Ak,2,function(z) { z/prod(z)^(1/length(z))})
		D = newD(D=D, Wk=Wk, Ak=Ak, G=G, d=d) 
		
		conv = c(testval(Wk=Wk,Ak=Ak,D=D,G=G), conv[1] )
		count = count +1
	}
#print(count)
	lam = numeric(G) 
	for (g in 1:G) lam[g] =sum(diag( D %*% diag(1/Ak[,g])%*% t(D) %*% Sk[,,g] ))/d
#print(apply(Ak,2,prod))
	
	val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
	for (g in 1:G) { 
		val$sigma[,,g]    = (D %*% diag(lam[g]*Ak[,g])%*% t(D))
		val$invSigma[,,g] = (D %*% diag(1/lam[g]*1/Ak[,g])%*% t(D))
		val$logdet[g]     =  d*log(lam[g])
		}
	return(val)
}



msVEE <- function(Sk=NULL, ng=NULL, eplison=1e-14, max.iter=100) {
	# Wk is a list of length G with matrix (p x p)
	# ng is a vector of length G of weights of Wk
	d = dim(Sk)[1];
	G = dim(Sk)[3];
	
	Wk = Sk
	W  = matrix(0,d,d)
	for (g in 1:G) {
		Wk[,,g] = Sk[,,g]*ng[g]
		W  = W + Wk[,,g]
	}
		
	C = W/det(W)^(1/d)
	invC = solve(C)
	lam = apply(Sk,3,function(z, invC) { sum(diag(z %*% invC)) }, invC=invC)/d

	val1 = sum(apply(Wk,3,function(z, invC) { sum(diag(z*invC)) }, invC= invC)/lam) +d*sum(ng*lam)
	conv = c(val1, Inf  )
	count = 1 
	while ( diff(conv)/conv[1]  > eplison & count < max.iter) {
#for (i in 1:max.iter) {
		C = sumSk.wt(Sk=Wk, wt=1/lam, d=d,G=G)
		C = C/det(C)^(1/d)		
		invC = solve(C)
		
		lam = apply(Sk,3,function(z, invC) { sum(diag(z %*% invC)) }, invC=invC)/d
		val1 = sum(apply(Wk,3,function(z, invC) { sum(diag(z*invC)) }, invC= invC )/lam) +d*sum(ng*lam)
		conv = c(val1, conv[1] )
		count = count +1
	}
#print(c(count,det(C), lam)	)
		
	invC = solve(C)
	val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
	for (g in 1:G) { 
		val$sigma[,,g]    = lam[g]*C
		val$invSigma[,,g] = 1/lam[g]*invC
		val$logdet[g]     = d*log(lam[g])
		}
	return(val)
}


msEVV <- function(Sk=NULL, ng=NULL, eplison=1e-12, max.iter= 100) {
	# Wk is a list of length G with matrix (p x p)
	# ng is a vector of length G of weights of Wk
	d = dim(Sk)[1];
	G = dim(Sk)[3];
	
	Wk = Sk
	Ck = Sk
	lam = numeric(G)
	for (g in 1:G) {
		Wk[,,g] = Sk[,,g]*ng[g]
		lam[g]    = det(Wk[,,g])^(1/d)
	}
	Ck = sweep(Wk, 3, 1/lam, FUN="*")	
	lam = sum(lam)
	
	
	val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
	for (g in 1:G) { 
		val$sigma[,,g]    = lam * Ck[,,g]
		val$invSigma[,,g] = 1/lam * solve(Ck[,,g])
		val$logdet[g]     =  d* log(lam ) 
		}
	return(val)
}





#        cat("N=","\n")
#        logl<-c(1:nmax)
#        cat(nrow(data),"\n")
