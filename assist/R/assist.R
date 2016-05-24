
###           ASSIST          
###
### Copyright 2004       Chunlei Ke <chunlei_ke@yahoo.com>,
###                      Yuedong Wang <yuedong@pstat.ucsb.edu>

alogit<- function(x) 1/(1+exp(-x))
addCorrGroup<- function(cor.struct, cor.form){
   pos.split<- regexpr("\\(", cor.struct)
   tmp.struct<-  substring(cor.struct, first=pos.split+1)
   pos.form<- regexpr("form=~", tmp.struct)
   if(pos.form==-1) pos.form<-1
   if(pos.form==1) {
        newStruct<- paste(substring(cor.struct, first=1, last=pos.split), 
        "form=~", cor.form, sep="")
       }
   else{
       newStruct<- substring(cor.struct, first=1, last=pos.split+pos.form-1)
       newStruct<- paste(newStruct,"form=~", cor.form, sep="")
    }      

   pos.formEnd<- regexpr(",", tmp.struct)
   if(pos.formEnd==-1) pos.formEnd<- regexpr(")", tmp.struct)

   newStruct<- paste(newStruct, substring(tmp.struct, first=pos.formEnd), 
               sep="")

   newStruct
}

bdiag<-function(x)
{
# construct block-diagonal matrix from a list of matrices
# Based on "bdiag" from Scott Chasalow 
#(scott.chasalow@cereon.com)
     if(is.matrix(x)) return(x)
     if(!is.list(x)) stop("dismatched input")
     n <- length(x)
     x <- lapply(x, function(y)
       if(length(y)) as.matrix(y) else stop("Zero-length component in x"))
     d <- array(unlist(lapply(x, dim)), c(2, n))
     rr <- d[1,  ]
     cc <- d[2,  ]
     rsum <- sum(rr)
     csum <- sum(cc)
     out <- array(0, c(rsum, csum))
     ind <- array(0, c(4, n))
     rcum <- cumsum(rr)
     ccum <- cumsum(cc)
     ind[1, -1] <- rcum[ - n]
     ind[2,  ] <- rcum
     ind[3, -1] <- ccum[ - n]
     ind[4,  ] <- ccum
     imat <- array(1:(rsum * csum), c(rsum, csum))
     iuse <- apply(ind, 2, function(y, imat)
     {
          imat[(y[1] + 1):y[2], (y[3] + 1):y[4]]
     }
     , imat = imat)
     iuse <- as.vector(unlist(iuse))
     out[iuse] <- unlist(x)
     out
}
chol.new<-
function(Q)
{
## to calculate decomposition Q=AA^T 
## Q is a symmetric matix
	if(max(abs(Q - t(Q))) < 0.0001) Q <- (Q + t(Q))/2 else stop("Q needs to be symmetric!")
	tmp <- eigen(Q)
	num <- sum(tmp$values < 1e-010)
	n <- ncol(Q) - num
	matDiag.prod(as.matrix(tmp$vector[, 1:n]), sqrt(tmp$values[1:n]), FALSE)
}
cubic<- 
function(s, t = s) 
{ 
## this function is to calculate the reproducing 
## kernel for cubic spline on [0,1] 
	n <- length(s) 
	m <- length(t) 
	if(any(c(s, t) > 1) || any(c(s, t) < 0)) 
		stop("cubic is only defined on [0, 1]") 
	matrix(.C("cubic_ker1", 
		as.double(s), 
		as.double(t), 
		as.integer(n), 
		as.integer(m), 
		val = double(n * m),
		PACKAGE="assist")$val, ncol = m, byrow = TRUE) 
} 
cubic2<- 
function(s, t = s) 
{ 
## reproducing  kernel for cubic spline on [0,T] 
	n <- length(s) 
	m <- length(t) 
        if(any(c(s, t) < 0)) 
                stop("cubic is only defined on [0, T]") 
	matrix(.C("cubic_ker2", 
		as.double(s), 
		as.double(t), 
		as.integer(n), 
		as.integer(m), 
		val = double(n * m),
		PACKAGE="assist")$val, ncol = m, byrow = TRUE) 
} 
dcrdr<-function(rkpk.obj, r)
{
## Splus inteface of dcrdr used for Bayesian CIs;
## rkpk.obj is an output from either disdr or gdsidr.
	ngrid <- ncol(r)
	if(is.null(rkpk.obj$swk))
		s <- rkpk.obj$s
	else s <- rkpk.obj$swk
	qraux <- rkpk.obj$qraux
	jpvt <- rkpk.obj$jpvt
	nobs <- rkpk.obj$nobs
	nNull <- rkpk.obj$nnull
        nnull<- max(nNull,1)
	if(is.null(rkpk.obj$qwk))
		q <- rkpk.obj$q
	else q <- rkpk.obj$qwk
	nlaht <- rkpk.obj$nlaht
	varht <- rkpk.obj$varht
	b <- varht/(10^nlaht)			
	storage.mode(s) <- "double"
	storage.mode(qraux) <- "double"
	storage.mode(jpvt) <- "integer"
	storage.mode(q) <- "double"
	storage.mode(nlaht) <- "double"
	storage.mode(varht) <- "double"
	storage.mode(r) <- "double"
	dcrdr <- .C("dcrdr_",
		s = s,
		lds = as.integer(nobs),
		nobs = as.integer(nobs),
		nnull = as.integer(nNull),
		qraux = qraux,
		jpvt = jpvt,
		q = q,
		ldq = as.integer(nobs),
		nlaht = nlaht,
		r = r,
		ldr = as.integer(nobs),
		nr = as.integer(ngrid),
		cr = double(nobs * ngrid),
		ldcr = as.integer(nobs),
		dr = double(nnull * ngrid),
		lddr = as.integer(nnull),
		wk = double(2 * nobs),
		info = integer(1),
		PACKAGE="assist")
	dcrdr
}


diagComp<-function(mat, vec)
{
## construct block diagonal matrix from
## a 'long' matrix and dim's of blocks 
   if(length(vec)==1) return(mat)
   len <- length(vec)
   cvec<-c(0, cumsum(vec))
   mlist<- as.list(1:len)
   mlist<- lapply(mlist, function(x, m, v) 
        if(nrow(m)>1) m[, (v[x]+1):v[x+1]] else matrix(m[, (v[x]+1):v[x+1]], nrow=1), 
        m=as.matrix(mat), v=cvec)
   bdiag(mlist)
}
dmudr<- 
function(y, q, s, weight = NULL, vmu = "v", theta = NULL, varht = NULL, tol = 0, 
	init = 0, prec = 1e-006, maxit = 30) 
{ 
## R interface of dmudr  
##  'weight', see dsidr 
## allow s=NULL	 
## handle options 
	if(is.null(vmu)) vmu <- "v" 
	if(is.null(init)) 
		init <- 0 
	if(is.null(tol)) 
		tol <- 0 
	if(is.null(prec)) 
		prec <- 1e-006 
	if(is.null(maxit)) 
		maxit <- 30 
	nobs <- length(y) 
	if(is.list(q)) 
		q <- array(unlist(q), dim = c(rep(length(y), 2), length(q))) 
	nq <- dim(q)[3] 
	if(!is.null(weight)) { 
		if(!is.null(s)) 
			s <- weight %*% s 
		for(i in 1:nq) 
			q[,  , i] <- weight %*% q[,  , i] %*% t(weight) 
		y1 <- weight %*% y 
	} 
	else y1 <- y 
	if(is.null(s)) { 
		s <- matrix(rep(0, nobs), ncol = 1) 
		nNull <- 0 
		nnull <- 1 
	} 
	else nNull <- nnull <- ncol(s) 
	if((init == 1) && (missing(theta) || is.null(theta))) 
		stop("initial values for theta is needed") 
	if(missing(theta) || is.null(theta)) 
		theta <- double(nq) 
	else theta <- as.double(theta) 
	storage.mode(s) <- "double" 
	storage.mode(q) <- "double" 
	storage.mode(y1) <- "double" 
	if(missing(varht) || is.null(varht)) { 
		if(vmu == "u") 
			stop("need varht for UBR method") 
		varht <- double(2) 
	} 
	else { 
		varht <- as.double(c(varht, 0)) 
	} 
	vmu2 <- (0:2)[vmu == c("v", "m", "u")] 
	result <- .C("dmudr_", 
		vmu = as.integer(vmu2), 
		s = s, 
		lds = as.integer(nobs), 
		nobs = as.integer(nobs), 
		nnull = as.integer(nNull), 
		q = q, 
		ldqr = as.integer(nobs), 
		ldqc = as.integer(nobs), 
		nq = as.integer(nq), 
		y = y1, 
		tol = as.double(tol), 
		init = as.integer(init), 
		prec = as.double(prec), 
		maxit = as.integer(maxit), 
		theta = theta, 
		nlaht = double(1), 
		score = double(1), 
		varht = as.double(varht), 
		c = double(nobs), 
		d = double(nnull), 
		wk = double(nobs * nobs * (nq + 2)), 
		info = integer(1),
		PACKAGE="assist") 
	if(result$info == -1) 
		stop("dimension error") 
	else if(result$info == -2) 
		stop("F_{2}^{T} Q_{*}^{theta} F_{2} !>= 0 ") 
	else if(result$info == -3) 
		stop("tuning parameters out of scope") 
	else if(result$info == -4) 
		stop("fails to converge within maxite steps") 
	else if(result$info == -5) 
		stop("fails to find a reasonable descent direction") 
	else if(result$info > 0) 
		stop("the matrix S is rank deficient: rank(S)+1") 
	result$fit <- s %*% result$d 
	tmp.result <- 0 
	for(i in 1:result$nq) 
		tmp.result <- tmp.result + (10^(result$theta[i]) * q[,  , i]) 
	if(vmu == "m") 
		result$varht <- result$varht[2] 
	else result$varht <- result$varht[1] 
	result$penalty <- sum(matVec.prod(tmp.result, result$c) * result$c) 
	result$fit <- result$fit + tmp.result %*% result$c 
	resultRss <- sum((y1 - result$fit)^2) 
	if(vmu=="u") result$df<- (length(y1)*result$score-resultRss)/(2*result$varht)
	else result$df <- length(y1) - resultRss/result$varht 
	
	if(!is.null(weight)) { 
		result$c <- t(weight) %*% result$c 
		result$fit <- solve(weight) %*% result$fit 
	} 
	result$resi <- y - result$fit 
	result$vmu <- vmu 
	result 
} 
dsidr<- 
function(y, q, s=NULL, weight = NULL, vmu = "v", varht = NULL, limnla = c(-10, 3),  
	job = -1, tol = 0) 
{ 
## R interface of dsidr with arguments as in the Fortran program  
## except that 'weight' is a matrix satisfying:   
##            (var(y))^(-1)=\sigma^2 * t(weight) %*% weight 
## s could be null	 
## handle options 
	if(is.null(vmu)) vmu <- "v" 
	if(is.null(job)) 
		job <- -1 
	if(is.null(tol)) 
		tol <- 0 
	if(is.null(limnla)) 
		limnla <- c(-10, 0) 
	if(!is.null(weight)) { 
		yy <- matVec.prod(weight, y, FALSE) 
		q <- weight %*% q %*% t(weight) 
		if(!is.null(s)) 
			s <- weight %*% s 
	} 
	else yy <- y 
	if(length(limnla) == 1) 
		limnla <- c(limnla, limnla) 
	if(missing(varht) || is.null(varht)) { 
		if(vmu == "u") 
			stop("need varht for UBR method") 
		varht <- double(2) 
	} 
	else { 
		varht <- as.double(c(varht, 0)) 
	} 
	nobs <- length(y) 
	if(is.null(s)) { 
		s <- matrix(rep(0, nobs), ncol = 1) 
		nNull <- 0 
		nnull <- 1 
	} 
	else nNull <- nnull <- ncol(s) 
	storage.mode(s) <- "double" 
	storage.mode(q) <- "double" 
	vmu2 <- (0:2)[vmu == c("v", "m", "u")] 
	result <- .C("dsidr_", 
		vmu = as.integer(vmu2), 
		s = s, 
		lds = as.integer(nobs), 
		nobs = as.integer(nobs), 
		nnull = as.integer(nNull), 
		y = as.double(yy), 
		q = q, 
		ldq = as.integer(nobs), 
		tol = as.double(tol), 
		job = as.integer(job), 
		limnla = as.double(limnla), 
		nlaht = double(1), 
		score = double(job + 2), 
		varht = varht, 
		c = double(nobs), 
		d = double(nnull), 
		qraux = double(nnull), 
		jpvt = integer(nnull), 
		wk = double(3 * nobs), 
		info = integer(1),
		PACKAGE="assist") 
	if(result$info == -1) 
		stop("dimension error") 
	else if(result$info == -2) 
		stop("F_{2}^{T} Q F_{2} !>= 0 ") 
	else if(result$info == -3) 
		stop("vmu is out of scope") 
	else if(result$info > 0) 
		stop("matrix S is rank deficient: rank(S)+1") 
      if(vmu == "m") 
		result$varht <- result$varht[2] 
	else result$varht <- result$varht[1] 
	result$penalty <- sum(matVec.prod(q, result$c) * result$c) 
	result$resi <- 10^(result$nlaht) * result$c 
	result$fit <- yy - result$resi 
	resultRss <- sum(result$resi^2) 
	if(result$job > 0){ 
		result$score <- result$score[ - length(result$score)] 
	}
      if(vmu=="u") result$df<- (nobs*min(result$score)-resultRss)/(2*result$varht)
	else result$df <- nobs - resultRss/result$varht
	if(!is.null(weight)) { 
		result$c <- t(weight) %*% result$c 
		result$fit <- solve(weight) %*% result$fit 
		result$resi <- y - result$fit 
	} 
	result$vmu <- vmu 
	result$nq <- 1 
	result 
} 
dsms<-function(rkpk.obj)
{
## Splus interface of dsms aiming to calculation of Bayesian CIs.
## rkpk.obj is a dsidr or gdsidr fit.
	if(is.null(rkpk.obj$swk)) ss <- rkpk.obj$s else ss <- rkpk.obj$swk
	qraux <- rkpk.obj$qraux
	jpvt <- rkpk.obj$jpvt
	nobs <- rkpk.obj$nobs
	nnull <- rkpk.obj$nnull
	if(is.null(rkpk.obj$qwk))
		q <- rkpk.obj$q
	else q <- rkpk.obj$qwk
	nlaht <- rkpk.obj$nlaht
	varht <- rkpk.obj$varht		      
	as.integer(nobs)
	as.integer(nnull)
	storage.mode(ss) <- "double"
	storage.mode(qraux) <- "double"
	storage.mode(jpvt) <- "integer"
	storage.mode(q) <- "double"
	storage.mode(nlaht) <- "double"
	storage.mode(varht) <- "double"	
	dsms <- .C("dsms_",
		s = ss,
		lds = nobs,
		nobs = nobs,
		nnull = nnull,
		jpvt = jpvt,
		q = q,
		ldq = nobs,
		nlaht = nlaht,
		sms = double(nnull * nnull),
		ldsms = nnull,
		wk = double(2 * nobs),
		info = integer(1),
		PACKAGE="assist")
	dsms
}

exponen<- 
function(s, t = s, r = 1) 
{ 
	n <- length(s) 
	m <- length(t) 
	matrix(.C("expLspline_ker", 
		as.double(s), 
		as.double(t), 
		as.double(r), 
		as.integer(n), 
		as.integer(m), 
		val = double(n * m),
		PACKAGE="assist")$val, ncol = m, byrow = TRUE) 
} 

gdmudr<- 
function(y, q, s, family, vmu = "v", varht = NULL, init = 0,  
	theta = NULL, tol1 = 0, tol2 = 0, prec1 = 1e-006, maxit1 = 30, prec2 =  
	1e-006, maxit2 = 30) 
{ 
## R inteface of dbmdr, dbimdr, dpmdr, dgmdr 
## weight is null.  
## handle options 
	if(is.null(vmu)) vmu <- "v" 
	if(is.null(tol1)) 
		tol1 <- 0 
	if(is.null(tol2)) 
		tol2 <- 0 
	if(is.null(prec1)) 
		prec1 <- 1e-006 
	if(is.null(prec2)) 
		prec2 <- 1e-006 
	if(is.null(maxit1)) 
		maxit1 <- 30 
	if(is.null(maxit2)) 
		maxit2 <- 30 
	if(family == "binomial") 
		nobs <- nrow(y) 
	else nobs <- length(y) 
	if(is.null(s)) { 
		s <- matrix(rep(0, nobs), ncol = 1) 
		nNull <- 0 
		nnull <- 1 
	} 
	else nNull <- nnull <- ncol(s) 
	if(is.list(q)) 
		q <- array(unlist(q), dim = c(rep(length(y), 2), length(q))) 
	nq <- dim(q)[3] 
	if(family == "binomial") 
		y <- t(y) 
	if((init == 1) && (missing(theta) || is.null(theta))) 
		stop("initial values for theta is needed") 
	if(missing(theta) || is.null(theta)) 
		theta <- double(nq) 
	else theta <- as.double(theta) 
	if(missing(varht) || is.null(varht)) 
		varht <- c(1, 0) 
	else varht <- c(varht, 0) 
	vmu2 <- (0:3)[vmu == c("v", "m", "u", "u~")] 
	storage.mode(s) <- "double" 
	storage.mode(q) <- "double" 
	storage.mode(y) <- "double" 
	switch(family, 
		binary = loadWhich <- "dbmdr_", 
		binomial = loadWhich <- "dbimdr_", 
		poisson = loadWhich <- "dpmdr_", 
		gamma = loadWhich <- "dgmdr_", 
		stop("unknown family")) 
	result <- .C(loadWhich, 
		vmu = as.integer(vmu2), 
		s = s, 
		lds = as.integer(nobs), 
		nobs = as.integer(nobs), 
		nnull = as.integer(nNull), 
		q = q, 
		ldqr = as.integer(nobs), 
		ldqc = as.integer(nobs), 
		nq = as.integer(nq), 
		y = y, 
		tol1 = as.double(tol1), 
		tol2 = as.double(tol2), 
		init = as.integer(init), 
		prec1 = as.double(prec1), 
		maxit1 = as.integer(maxit1), 
		prec2 = as.double(prec2), 
		maxit2 = as.integer(maxit2), 
		theta = theta, 
		nlaht = double(1), 
		score = double(1), 
		varht = as.double(varht), 
		c = double(nobs), 
		d = double(nnull), 
		eta = double(nobs), 
		wk = double(nobs * nobs * (nq + 2)), 
		swk = double(nnull * nobs), 
		qwk = double(nobs * nobs * nq), 
		ywk = double(nobs), 
		u = double(nobs), 
		w = double(nobs), 
		info = integer(1),
		PACKAGE="assist") 
	if(result$info > 0) 
		stop("the matrix S is rank deficient: rank(S)+1")
	else { 
		if(result$info < 0) 
			switch( - result$info, 
				stop("dimension error"), 
				stop("F_{2}^{T} Q_{*}^{theta} F_{2} !>= 0"), 
				stop("tuning parameters are out of scope"), 
				stop( 
				  "DMUDR fails to converge within maxiter1 steps" 
				  ), 
				stop( 
				  "fails to find a reasonable descent direction" 
				  ), 
				stop( 
				  "Logit estimation fail to converge within maxiter2 steps" 
				  ), 
				stop("There are some w's equals to zero")) 
	} 
	if(vmu == "m") 
		result$varht <- result$varht[2] 
	else result$varht <- result$varht[1] 
	result$rss <- (10^result$nlaht * result$c)/sqrt(result$w) 
	result$rss <- sum(result$rss^2) 
	if(vmu=="u") result$df<- (nobs*result$score-result$rss)/(2*result$varht)
	else result$df <- nobs - result$rss/result$varht 
	
	switch(family, 
		binary = result$fit <- alogit(result$eta), 
		binomial = { 
			result$fit <- alogit(result$eta) 
			result$fit <- y[1,  ] * result$fit 
		} 
		, 
		poisson = result$fit <- exp(result$eta), 
		gamma = result$fit <- exp(result$eta)) 
	result$vmu <- vmu 
	result 
} 
gdsidr<- 
function(y, q, s, family, vmu = "v", varht = NULL, limnla = c( 
	-10, 3), maxit = 30, job = -1, tol1 = 0, tol2 = 0, prec = 1e-006) 
{ 
##  R interface of dbsdr, dbisdr, dpsdr and dgsdr; 
##  weight is NULL. 
##  allow S=NULL  
## handle options 
	if(is.null(vmu)) vmu <- "v" 
	if(is.null(job)) 
		job <- -1 
	if(is.null(tol1)) 
		tol1 <- 0 
	if(is.null(tol2)) 
		tol2 <- 0 
	if(is.null(prec)) 
		prec <- 1e-006 
	if(is.null(maxit)) 
		maxit <- 30 
	if(is.null(limnla)) 
		limnla <- c(-10, 3) 
	if(family == "binomial") { 
		nobs <- nrow(y) 
		y <- t(y) 
	} 
	else nobs <- length(y) 
	if(is.null(s)) { 
		s <- matrix(rep(0, nobs), ncol = 1) 
		nNull <- 0 
		nnull <- 1 
	} 
	else nNull <- nnull <- ncol(s) 
	storage.mode(s) <- "double" 
	storage.mode(q) <- "double" 
	storage.mode(y) <- "double" 
	switch(family, 
		binary = loadWhich <- "dbsdr_", 
		binomial = loadWhich <- "dbisdr_", 
		poisson = loadWhich <- "dpsdr_", 
		gamma = loadWhich <- "dgsdr_", 
		stop("unknown family")) 
	if(length(limnla) == 1) 
		limnla <- c(limnla, limnla) 
	if(missing(varht) || is.null(varht)) 
		varht <- c(1, 0) 
	else varht <- c(varht, 0) 
	vmu2 <- (0:3)[vmu == c("v", "m", "u", "u~")] 
	result <- .C(loadWhich, 
		vmu = as.integer(vmu2), 
		s = s, 
		lds = as.integer(nobs), 
		nobs = as.integer(nobs), 
		nnull = as.integer(nNull), 
		y = y, 
		q = q, 
		ldq = nobs, 
		tol1 = as.double(tol1), 
		tol2 = as.double(tol2), 
		job = as.integer(job), 
		limnla = as.double(limnla), 
		prec = as.double(prec), 
		maxit = as.integer(maxit), 
		nlaht = double(1), 
		score = double(job + 2), 
		varht = as.double(varht), 
		c = double(nobs), 
		d = double(nnull), 
		eta = double(nobs), 
		qraux = double(nnull), 
		jpvt = integer(nnull), 
		wk = double(3 * nobs), 
		swk = double(nnull * nobs), 
		qwk = double(nobs * nobs), 
		ywk = double(nobs), 
		u = double(nobs), 
		w = double(nobs), 
		info = integer(1),
		PACKAGE="assist") 
	if(result$info > 0) 
		stop("matrix S is rank deficient: rank(S)+1") 
	else { 
		if(result$info < 0) 
			switch( - result$info, 
				stop("dimension error"), 
				stop("F_{2}^{T} Q F_{2} !>= 0"), 
				stop("vmu is out of scope"), 
				stop("fail to converge at the maxiter steps"), 
				warning("There are some w's equals to zero")) 
	} 
	switch(family, 
		binary = result$fit <- alogit(result$eta), 
		binomial = { 
			result$fit <- alogit(result$eta) 
			result$fit <- y[1,  ] * result$fit 
		} 
		, 
		poisson = result$fit <- exp(result$eta), 
		gamma = result$fit <- exp(result$eta)) 
	if(result$job > 0) result$score <- result$score[ - length(result$score)]	 
## needed to add something above this line ## 
	if(vmu == "m") 
		result$varht <- result$varht[2] 
	else result$varht <- result$varht[1] 
	result$resi <- (10^result$nlaht * result$c)/sqrt(result$w) 
	resultRss <- sum(result$resi^2) 
	if(vmu=="u") result$df<- (nobs*min(result$score)-resultRss)/(2*result$varht)
	else result$df <- nobs - resultRss/result$varht 

	result$vmu <- vmu 
	result$nq <- 1 
	result 
} 
getFunInfo<- function(object){
## extract function information
## used in SMR, SNM and nnr
## allow general forms
  if(!is.list(object)) object<- list(object)
  lengthF<- length(object)
  fName<- NULL
  fNullForm<- fRkForm<- f.arg<-list()
  for(i in 1:lengthF){
    compFun<-  splitFormula(getResponseFormula(object[[i]]), sep="+")
    fName<- c(fName, unlist(lapply(compFun, function(x) as.character(x[[2]])[1])))
    f.arg<- c(f.arg, lapply(compFun, function(x) as.character(x[[2]])[-1]))
    f.comp<- object[[i]][[3]]
    if(length(f.comp)==2){ 
        fRkForm<- c(fRkForm, rep(list(f.comp[[2]]), length(compFun)))        
       }
    else if(is.null(f.comp$nb) || is.null(f.comp$rk)){
        fRkForm<- c(fRkForm, rep(list(f.comp[[3]]), length(compFun)))
        for(IndexNull in 1:length(compFun)) 
               fNullForm[[length(f.arg)-IndexNull+1]]<-f.comp[[2]]
      }
   }
  list(fName=fName, fNullForm=fNullForm, fRkForm=fRkForm, f.arg=f.arg)
}
getParaModelMatrix<-function(params, data, nobs)
{
## extract parameter names and model.matrix
## used internally
   matrix.para<- NULL
   length.para<- NULL
   para.name<- NULL
   if(is.call(params)) params<- eval(params)
   if(!is.list(params) ) params<- list(params)
   for(i in 1:length(params)){         
     if(as.character(getCovariateFormula(params[[i]])[[2]])[[1]]=="1"){
        temp.matrix<-matrix(rep(1, nobs), ncol=1)
       }
     else{
         temp.matrix<- model.matrix2(getCovariateFormula(params[[i]]), data)
         }
     temp.fixed<-splitFormula(getResponseFormula(params[[i]]), sep="+")
     temp.fixed<- unlist(lapply(temp.fixed, function(x) as.character(x[[2]])))

     para.name<-c(para.name, temp.fixed )
     length.para<- c(length.para, rep(ncol(temp.matrix),length(temp.fixed)))

     if(nrow(temp.matrix)==1) 
        temp.matrix<- matrix(rep(temp.matrix, nobs), nrow=nobs, byrow=TRUE)

     temp.matrix<-matrix(rep(temp.matrix, length(temp.fixed)),
         nrow=nobs, byrow=FALSE)
     matrix.para<- cbind(matrix.para, temp.matrix)
   }
   list(para.name=para.name, matrix.para=matrix.para, length.para=length.para)
}

getParaValue<- function(length.fix, length.random, matrix.fix, matrix.random, 
               start.fixed, start.random, para, para.random, nobs, para.fixed, nobs.in)
{
## used for snm
## not for general purpose  
  startFixMat<- diagComp(matrix(start.fixed, nrow=1),length.fix)
  para.v<- matrix.fix%*%t(startFixMat)
 
  if((length(para)-length(para.fixed))>0)
       para.v <- cbind(para.v, 
            matrix( rep(0, (length(para)-length(para.fixed))*nobs), nrow=nrow(para.v)))

  if(!is.list(nobs.in)) nobs.in<-list(nobs.in)

  sumRand<- 0
  for(i in 1:length(start.random)){  
     tmpNobs<- nobs.in[[i]]     
     sumRand<-sumRand+ apply(start.random[[i]], 2,
           function(x, nobs.in) rep(x, nobs.in), nobs.in=tmpNobs) 
#          function(x, nobs.in) rep(x, nobs.in), nobs.in=nobs.in[[i]])
   }

 start.random<- matrix.random[, 1:ncol(sumRand)]*sumRand
 lengthRan<- rep(1:length(length.random), length.random)
 para.v.random<-t(rowsum(t(start.random), lengthRan))
 
 repVar<- match(para.random, para)
 para.v[, repVar]<-para.v[,repVar]+para.v.random    

 para.v<- data.frame(para.v)
 names(para.v)<- para
  para.v
}


hat.ssr<-
function(ssr.obj)
{
## To calculate jat matrix for ssr.object ##
	if(length(ssr.obj$q) > 1) is.dmudr <- TRUE else is.dmudr <- FALSE
	if(is.dmudr)
		theta <- 10^(ssr.obj$rkpk.obj$theta)
	else theta <- 1
	nobs <- nrow(ssr.obj$s)	
# call (g)dsidr fit for dmudr
	if(is.dmudr) {
		qSum <- 0
		for(i in 1:length(theta)) {
			qSum <- qSum + ssr.obj$q[[i]] * theta[i]
		}
		swk <- ssr.obj$s
		y <- ssr.obj$y
		if(is.null(ssr.obj$expand.call$family) || ssr.obj$expand.call$
			family == "gaussian") {
			dsidr.fit <- dsidr(swk, q = qSum, y = y, weight = 
				ssr.obj$weight, vmu = ssr.obj$call$spar, tol
				 = ssr.obj$call$control$tol, job = ssr.obj$call$
				control$job, limnla = ssr.obj$call$limnla)
		}
		else {
			dsidr.fit <- gdsidr(s = swk, q = qSum, y = y, family = 
				ssr.obj$expand.call$family, vmu = ssr.obj$
				expand.call$spar, tol1 = ssr.obj$expand.call$
				control$tol, tol2 = ssr.obj$expand.call$control$
				tol.g, maxit = ssr.obj$expand.call$control$
				maxit.g, job = ssr.obj$expand.call$control$job, 
				limnla = ssr.obj$expand.call$limnla, prec = 
				ssr.obj$expand.call$control$prec)
		}
	}
	else dsidr.fit <- ssr.obj$rkpk.obj
	cr <- 10^(dsidr.fit$nlaht) * matrix(dcrdr(dsidr.fit, diag(nobs))$cr, 
		ncol = nobs, byrow = FALSE)
	diag(cr) <- diag(cr - 1)
 - cr
}

ident<-
function(x, y = x)
{
## scale y by: (y-min(x))/(max(x)-min(x)) ##
## this is used in ssr                    ##
	if(is.vector(x) && !is.list(x)) {
			if(!is.vector(y) || is.list(y))
				stop(" x and y must be of the same type")
			else{                          
                          if(is.numeric(resul<- y)) resul <- (y - min(x))/(max(x) - min(x))
                          else warning("Non-numeric was not scaled")
                         }
	}	
	else if(is.data.frame(x)) {
                        if(!is.data.frame(y)) stop(" x and y doesn't match")
			resul <- NULL
			for(nn in names(y)) {
                                if(is.numeric(y[,nn])){
				    resul <- cbind(resul, (y[, nn] - min(x[, nn]))/(
				  max(x[, nn]) - min(x[, nn])))
                                 }
                                else{
                                   resul<- cbind(resul, y[, nn])
                                   warning("Non-numeric variable was not scaled")
                                  }
			}
			resul <- data.frame(resul)
			names(resul) <- names(y)
		}
	else if(is.matrix(x)) { 		
			if(!is.matrix(y))
		           stop("x and y must be of the same type")
			else {
                            if(is.numeric(resul<- y)){
			    resul <- NULL
			    for(i in 1:ncol(x))
				      resul <- cbind(resul, (y[, i] - min(x[, i
				        ]))/(max(x[, i]) - min(x[, i])))
                            }
                            else warning("Non-numeric variable was not scaled")
			}
		}			
	else if(is.list(x)) {
		        if(!is.list(y) || (length(y) != length(x)))
			    stop(" x and y must be of the same type")
		        else {
				    resul <- list(length(y))
				    for(i in 1:length(x))
                                      if(is.numeric(resul[[i]]<-y[[i]]))
				          resul[[i]] <- (y[[i]] - min(x[[i]]))/(max(
				                x[[i]]) - min(x[[i]]))
                                      else warning("non-numeric variable was not scaled")
		         }
	}
        else if(!is.numeric(resul<-y))           
           warning("non-numeric variable was not scaled")           
        else stop("unknown input type")						
	resul
}

kron<- 
function(x, y = x) 
{ 
## this function is to construct kernels for basis functions  ## 
## for example kron(list(f1, f2))=f1*f1^T+f2*f2^T             ## 
	if(!is.list(x)) matrix(kronecker(x, y), ncol = length(y), byrow=TRUE) else { 
		resul <- 0 
		for(i in 1:length(x)) { 
			if(length(x[[i]]) == 1) 
				resul <- resul + x[[i]] * y[[i]] 
			else resul <- resul + matrix(kronecker(x[[i]], y[[i]]), byrow=TRUE, 
				  ncol = length(y[[i]])) 
		} 
		resul 
	} 
} 
linSinCos<- 
function(s, t = s) 
{ 
	n <- length(s) 
	m <- length(t) 
	matrix(.C("linPeriod_ker", 
		as.double(s), 
		as.double(t), 
		as.integer(n), 
		as.integer(m), 
		val = double(n * m),
		PACKAGE="assist")$val, ncol = m, byrow = TRUE) 
} 
linear<-
function(s, t = s)
{
## reproducing kernel for cubic spline on [0,1]
        n <- length(s)
        m <- length(t)
        if(any(c(s, t) > 1) || any(c(s, t) < 0))
                stop("linear is only defined on [0, 1]")
        matrix(.C("linear_ker1",
                as.double(s),
                as.double(t),
                as.integer(n),
                as.integer(m),
                val = double(n * m),
		PACKAGE="assist")$val, ncol = m, byrow = TRUE)
}
linear2<-
function(s, t = s)
{
        n <- length(s)
        m <- length(t)
        if(any(c(s, t) < 0))
                stop("linear is only defined on [0, T]")
        matrix(.C("linear_ker2",
                as.double(s),
                as.double(t),
                as.integer(n),
                as.integer(m),
                val = double(n * m),
		PACKAGE="assist")$val, ncol = m, byrow = TRUE)
}
list.prod<- 
function(x, y, fun="*") 
{ 
## operation of two lists/vector(s)     
## of the same length                                
   if(length(x) != length(y)) stop("Unequal length!") 
   lapply(as.list(1:length(y)), function(x, x1, y1, fun) 
        do.call(fun, list(x1[[x]], y1[[x]])), x1=x, y1=y, fun=fun)  
} 
 
logitKer<-
function(s, t = s)
{
        n <- length(s)
        m <- length(t)
        matrix(.C("logit_ker",
                as.double(s),
                as.double(t),
                as.integer(n),
                as.integer(m),
                val = double(n * m),
		PACKAGE="assist")$val, ncol = m, byrow = TRUE)
}
lspline<-
function(x, y = x, type = "exp", ...)
{
	call<- match.call()
	type<- as.character(call$type)
	switch(type,
		exp = exponen(x, y, ...),
		logit = logitKer(x, y),
		sine0 = sine0(x, y),
		sine1 = sine1(x, y),
		linSinCos = linSinCos(x, y),
		stop("unknown input of type"))
}
matDiag.prod<- 
function(mat, vec, left = TRUE) 
{ 
## calculate diag(vec)%*%mat or mat%*%vec ##  
	if(left) { 
		if(length(vec) != nrow(mat)) 
			stop("length of Vector does not match row of Matrix") 
                mat*vec 
	} 
	else { 
		if(length(vec) != ncol(mat)) 
			stop("length of Vector does not match column of Matrix" 
				) 
                t(t(mat)*vec) 
	} 
} 
matVec.prod<- 
function(mat, vec, left = TRUE) 
{ 
## to multiply mat with vec 
## v%*%M or M%*%v  
	if(left && (length(vec) != nrow(mat))) stop( 
			"length of Vector does not match row of Matrix") 
	if(!left && (length(vec) != ncol(mat))) 
		stop("length of Vector does not match column of Matrix") 
	as.vector(apply(mat, 1 + left, function(x, z) 
	sum(x * z), z = vec)) 
} 
model.matrix2<- 
function(formula, data = list(), ...) 
{ 
## enhanced model.matrix to be used in snr, snm, nnr  
   coForm <- getCovariateFormula(formula) 
   if((length(coForm[[2]])==1)&& (as.character(coForm[[2]]) == "1") && !missing(data)) { 
	val <- as.matrix(data.frame(rep(1, nrow(data)))) 
	dimnames(val)[[2]] <- "(Intercept)" 
	val 
   } 
  else model.matrix(eval(formula), data = data, ...) 
} 

paste2<- function(nam, value){ 
## return a string of the form "a=1, b=2" 
## used in rkEval 
    x<- paste(nam, "=", value, sep="") 
    x[nam==""]<-paste(value)[nam==""] 
    re<- x[1] 
    if(length(nam)>1) { 
       for(len in 2:length(nam)) re<- paste(re, x[len], sep=",") 
     } 
   re 
} 
periodic<- 
function(s, t = s, order = 2) 
{ 
        s <- s %% 1 
        t <- t %% 1 
        n <- length(s) 
        m <- length(t) 
        if(order == 2) { 
                matrix(.C("period_ker", 
                        as.double(s), 
                        as.double(t), 
                        as.integer(n), 
                        as.integer(m), 
                        val = double(n * m),
			PACKAGE="assist")$val, ncol = m, byrow = TRUE) 
        } 
        else { 
                matrix(.C("mperiod_ker", 
                        as.double(s), 
                        as.double(t), 
                        as.integer(n), 
                        as.integer(m), 
                        as.integer(order), 
                        val = double(n * m),
			PACKAGE="assist")$val, ncol = m, byrow = TRUE) 
        } 
} 

quintic<-
function(s, t = s)
{
	n <- length(s)
	m <- length(t)
	if(any(c(s, t) > 1) || any(c(s, t) < 0))
		stop("quintic is only defined on [0, 1]")
	matrix(.C("quintic_ker1",
		as.double(s),
		as.double(t),
		as.integer(n),
		as.integer(m),
		val = double(n * m),
		PACKAGE="assist")$val, ncol = m, byrow = TRUE)
}
quintic2<-
function(s, t = s)
{
        n <- length(s)
        m <- length(t)
        if(any(c(s, t) < 0))
                stop("quintic is only defined on [0, T]")
        matrix(.C("quintic_ker2",
                as.double(s),
                as.double(t),
                as.integer(n),
                as.integer(m),
                val = double(n * m),
		PACKAGE="assist")$val, ncol = m, byrow = TRUE)
}
rk.prod<-
function(x, ...)
{
    rkProd<-function(s, t){
	if(is.matrix(s)) n <- ncol(s) else n <- length(s)
	if(is.matrix(t)) m <- ncol(t)
	else m <- length(t)
	if(n != m) stop("the dims of inputs must match")
	if(!is.matrix(s)) s <- kronecker(s, s)
	if(!is.matrix(t)) t <- kronecker(t, t)
	matrix(as.vector(s) * as.vector(t), ncol = n)
   }
   resul <- x
   for(y in list(...))
	resul <- rkProd(resul, y)
   resul
}
rkEval<-function(rk, g, g2) 
{ 
## eval rks 
## used internally  
   if(!is.call(rk)) rk<- as.call(rk)[[1]] 
   if(as.character(rk[[1]])!='list') rk<- list(NULL, rk) 
   leng <- length(rk) - 1 
   resul <- list(leng) 
   for(i in 1:leng) {    
      if(rk[[i + 1]][[1]] != "rk.prod") { 
         arg.tmp<- as.character(rk[[i+1]]) 
         arg.name<- names(rk[[i+1]]) 
         if(length(arg.tmp)< 3){          
	    resul[[i]] <- paste(arg.tmp[1], "(eval(expression(", arg.tmp[2],  
                    "), g),eval(expression(", arg.tmp[2],"), g2))", sep="") 
          } 
      else{             
        resul[[i]] <- paste(arg.tmp[1], "(eval(expression(", arg.tmp[2],  
                    "), g),eval(expression(", arg.tmp[2],"), g2),",sep="") 
        resul[[i]]<- paste(resul[[i]], paste2(arg.name[-c(1,2)], arg.tmp[-c(1, 2)]), ")", 
                 sep="") 
           }                        
      resul[[i]] <- eval(parse(text = resul[[i]]))	 
     } 
     else { 
	leng2 <- length(rk[[i + 1]]) 
        arg.tmp<- as.character(rk[[i+1]][[2]]) 
        arg.name<- names(rk[[i+1]][[2]]) 
        if(length(arg.tmp)< 3){          
	    resul[[i]] <- paste(arg.tmp[1], "(eval(expression(", arg.tmp[2],  
                    "), g),eval(expression(", arg.tmp[2],"), g2))", sep="") 
           } 
      else{ 
        resul[[i]] <- paste(arg.tmp[1], "(eval(expression(", arg.tmp[2],  
                    "), g),eval(expression(", arg.tmp[2],"), g2),",sep="") 
        resul[[i]]<- paste(resul[[i]], paste2(arg.name[-c(1,2)], arg.tmp[-c(1, 2)]), ")", 
                 sep="") 
            } 
	resul[[i]] <- eval(parse(text = resul[[i]])) 
	for(j in 3:leng2) { 
           arg.tmp<- as.character(rk[[i+1]][[j]]) 
           if(length(arg.tmp)< 3){          
	        temp1 <- paste(arg.tmp[1], "(eval(expression(", arg.tmp[2],  
                    "), g),eval(expression(", arg.tmp[2],"), g2))", sep="") 
           } 
           else{ 
               temp1<- paste(arg.tmp[1], "(eval(expression(", arg.tmp[2],  
                    "), g),eval(expression(", arg.tmp[2],"), g2),",sep="") 
               temp1<- paste(temp1, paste2(arg.name[-c(1,2)], arg.tmp[-c(1, 2)]), ")", 
                 sep="") 
 
          }	 
	temp1 <- eval(parse(text = temp1)) 
	resul[[i]] <- rk.prod(resul[[i]], temp1)				     
	} 
      } 
    } 
  resul 
} 
septic<-
function(s, t = s)
{
        n <- length(s)
        m <- length(t)
        if(any(c(s, t) > 1) || any(c(s, t) < 0))
                stop("septic is only defined on [0, 1]")
        matrix(.C("septic_ker1",
                as.double(s),
                as.double(t),
                as.integer(n),
                as.integer(m),
                val = double(n * m),
		PACKAGE="assist")$val, ncol = m, byrow = TRUE)
}
septic2<-
function(s, t = s)
{
        n <- length(s)
        m <- length(t)
        if(any(c(s, t) < 0))
                stop("septic is only defined on [0, T]")
        matrix(.C("septic_ker2",
                as.double(s),
                as.double(t),
                as.integer(n),
                as.integer(m),
                val = double(n * m),
		PACKAGE="assist")$val, ncol = m, byrow = TRUE)
}
shrink0<-
function(x, y = x)
{
	m <- length(y)
	n <- length(x)
	x <- unclass(as.factor(x))
	y <- unclass(as.factor(y))
	val <- .C("factor_ker",
		as.integer(x),
		as.integer(y),
		as.integer(n),
		as.integer(m),
		val = integer(n * m),
		PACKAGE="assist")$val
	matrix(val, ncol = m, byrow = TRUE)
}
shrink1<-
function(x, y = x)
{
	m <- length(y)
	n <- length(x)
	x <- unclass(as.factor(x))
	y <- unclass(as.factor(y))
	val <- .C("factor_ker",
		as.integer(x),
		as.integer(y),
		as.integer(n),
		as.integer(m),
		val = integer(n * m),
		PACKAGE="assist")$val
	matrix(val, ncol = m, byrow = TRUE) - 1/length(unique(x))
}
sine0<-
function(s, t = s)
{
	n <- length(s)
	m <- length(t)
	matrix(.C("sinLspline_ker0",
		as.double(s),
		as.double(t),
		as.integer(n),
		as.integer(m),
		val = double(n * m),
		PACKAGE="assist")$val, ncol = m, byrow = TRUE)
}
sine1<-
function(s, t = s)
{
	n <- length(s)
	m <- length(t)
	matrix(.C("sinLspline_ker1",
		as.double(s),
		as.double(t),
		as.integer(n),
		as.integer(m),
		val = double(n * m),
		PACKAGE="assist")$val, ncol = m, byrow = TRUE)
}

sine4p<-
function(s, t = s)
{
	n <- length(s)
	m <- length(t)
	matrix(.C("sinLspline_ker4p",
		as.double(s),
		as.double(t),
		as.integer(n),
		as.integer(m),
		val = double(n * m),
		PACKAGE="assist")$val, ncol = m, byrow = TRUE)
}
sphere<- 
function(x, y = x, order = 2) 
{ 
## calculate RK for pseudo spherical spline  
	if(order < 2 || order > 6) stop("order can only be 2, 3, 4, 5, 6") 
	if(is.matrix(x) && is.matrix(y)) { 
		x1 <- x[, 1] 
		x2 <- x[, 2] 
		y1 <- y[, 1] 
		y2 <- y[, 2] 
	} 
	else { 
		if(is.list(x) && is.list(y)) { 
			x1 <- x[[1]] 
			x2 <- x[[2]] 
			y1 <- y[[1]] 
			y2 <- y[[2]] 
		} 
	} 
	N <- length(x1) 
	M <- length(y1) 
	matrix(.C("sphere_ker", 
		as.double(x1), 
		as.double(x2), 
		as.double(y1), 
		as.double(y2), 
		as.integer(N), 
		as.integer(M), 
		as.integer(order), 
		val = double(N * M),
		PACKAGE="assist")$val, ncol = M, byrow = TRUE) 
} 

ssr.control<- 
function(job = -1, tol = 0, init = 0, theta = NULL, prec = 1e-006, maxit = 30,  
	tol.g = 0, prec.g = 1e-006, maxit.g = 30) 
{ 
	if((init == 1) && is.null(theta)) 
		stop("initial values for theta is needed") 
	list(job = job, tol = tol, init = init, theta = theta, prec = prec,  
		maxit = maxit, tol.g = tol.g, prec.g = prec.g, maxit.g =  
		maxit.g) 
} 
ssr<-
function (formula, rk, data = sys.parent(), subset, weights = NULL, 
    correlation = NULL, family = "gaussian", scale = FALSE, spar = "v", 
    varht = NULL, limnla = c(-10, 3), control = list()) 
{
    Call <- match.call()
    m <- match.call(expand.dots = FALSE)
    m$rk <- m$family <- m$correlation <- m$scale <- m$spar <- m$weights <- m$control <- m$scale <- m$varht <- m$limnla <- m$theta <- NULL
    m$na.action <- na.fail
    m[[1]] <- as.name("model.frame")
    m1 <- eval(m, sys.parent())
    Terms <- attr(m1, "terms")
    y <- model.extract(m1, "response")
#  y<- eval(getResponseFormula(formula)[[2]], data)
    if (family == "binomial") {
        y <- cbind(y[, 1] + y[, 2], y[, 1])
        n <- nrow(y)
    }
    else n <- length(y)
    varName <- unique(all.vars(Call$rk))
    dataOrg <- data
    if (eval(m$formula)[[3]]=="-1") S<- NULL 
    else if (scale) {
        if (missing(data)) {
            data <- as.list(varName)
            data <- data.frame(lapply(data, get))
            names(data) <- varName
        }
        for (var in intersect(varName, names(data))) data[[var]] <- ident(data[[var]])
        m$data <- data
        m <- eval(m, sys.parent())
        Terms <- attr(m, "terms")
        S <- model.matrix(Terms, m, NULL)
    }
    else S <- model.matrix(Terms, m1, NULL)
    if (length(S) == 0) 
        S <- NULL
    if ((missing(correlation) || is.null(correlation)) && (missing(weights) || 
        is.vector(weights <- eval(weights, data)) || is.matrix(weights))) {
        if (missing(weights) || is.null(weights)) 
            weights <- NULL
        else {
            if (!(is.vector(weights) || is.matrix(weights))) 
                stop("input of weight does not match!")
            if (is.vector(weights)) {
                if (length(weights) != n) 
                  stop("length of weight must match that of y")
                else weights <- diag(1/sqrt(weights))
            }
            else {
                if (is.matrix(weights) && any(dim(weights) != 
                  n)) 
                  stop("dim of weight must match with x and y")
                else weights <- solve(t(chol(weights)))
            }
        }
    }
    cor.est <- var.est <- NULL
    Q <- eval(Call$rk, data)
    if (is.matrix(Q)) {
        Q <- list(Q)
        Call$rk <- list(Call$rk)
    }
    if (!missing(subset)) {
        for (comp in 1:length(Q)) {
            Q[[comp]] <- Q[[comp]][eval(Call$subset, envir = data), 
                eval(Call$subset, envir = data)]
        }
    }
    ctrl.vals <- ssr.control()
    if (!missing(control)) {
        for (i in names(control)) ctrl.vals[[i]] <- control[[i]]
    }
    is.lme <- 0
    lme.obj <- NULL
    if (family == "gaussian") {
        if ((missing(correlation) || is.null(correlation)) && 
            (missing(weights) || is.vector(weights) || is.matrix(weights) || 
                is.null(weights))) {
            if (length(Q) > 1) {
                rk.obj <- dmudr(y = y, q = Q, s = S, weight = weights, 
                  vmu = spar, varht = varht, init = ctrl.vals$init, 
                  prec = ctrl.vals$prec, tol = ctrl.vals$tol, 
                  theta = ctrl.vals$theta, maxit = ctrl.vals$maxit)
##                lambda <- 10^rk.obj$theta
		lambda<- 10^(rk.obj$nlaht-rk.obj$theta)
            }
            else {
                rk.obj <- dsidr(y = y, q = Q[[1]], s = S, weight = weights, 
                  vmu = spar, varht = varht, limnla = limnla, 
                  job = ctrl.vals$job, tol = ctrl.vals$tol)
                lambda <- 10^(rk.obj$nlaht)
            }
        }
        else {
            if (is.null(S)) 
                stop("Empty null-space is not allowed!")
            tmpGrp<- rep(1, nrow(S))
            if (!is.null(correlation)) {
#print(class(correlation))
#                cor.form <- formula.corStruct(correlation)
		cor.form <- formula(correlation)
                if (!is.null(cor.form)) 
                  cor.group <- getGroupsFormula(cor.form)
                else cor.group <- NULL
                if (!is.null(cor.group)) {
                  cor.group <- as.character(cor.group)[[2]]
                  cor.group <- paste("tmpGrp", cor.group, sep = "/")
                  cor.form.new <- paste(as.character(getCovariateFormula(cor.form))[[2]], 
                    "|", cor.group, sep = "")
                  cor.form.new <- addCorrGroup(deparse(Call$cor), 
                    cor.form.new)
                  correlation <- eval(parse(text = cor.form.new))
                }
                else cor.form.new <- NULL
            }
            if (!(is.null(correlation) || is.null(cor.form.new))) {
                grpfac <- unique(all.vars(getGroupsFormula(cor.form)))
                if (!is.null(grpfac)) {
                  if (is.null(data)) {
                    dataCollect <- NULL
                    for (i in grpfac) dataCollect <- cbind(dataCollect, 
                      get(i))
                    dataCollect <- data.frame(dataCollect)
                    names(dataCollect) <- grpfac
                    grpOrder <- order.rows(dataCollect, grpfac)
                  }
                  else grpOrder <- order.rows(data, grpfac)
                }
                else grpOrder <- NULL
            }
            else grpOrder <- NULL
            if (!is.null(grpOrder)) {
                ord.tmp <- 1:length(grpOrder)
                grpOrder2 <- match(ord.tmp, grpOrder)
            }
            else grpOrder2 <- NULL
            if (spar != "m") 
                stop("Only GML method is available")
            tmp.z <- lapply(Q, chol.new)
            if (length(Q) < 2) {
                tmp.z <- tmp.z[[1]]
                if (!missing(data)) {
                  lmeData <- data
                  lmeData$S <- S
                  lmeData$tmp.z <- tmp.z
		  lmeData$tmpGrp<- tmpGrp
                  lme.obj <- lme(formula, random = list(tmpGrp=pdIdent(~tmp.z - 
                    1)), correlation = correlation, weights = weights, 
                    data = lmeData)
                }
                else {
                  lme.obj <- lme(formula, random = list(tmpGrp= pdIdent(~tmp.z - 
                    1)), correlation = correlation, weights = weights)
                }
            }
            else {
                tmp.block <- NULL
                tmp.num <- rep(0, n)
                tmp.zz <- NULL
                for (i in 1:length(tmp.z)) {
                  tmp.zz <- cbind(tmp.zz, tmp.z[[i]])
                  tmp.num[i + 1] <- ncol(tmp.z[[i]]) + tmp.num[i]
                  tmp <- list(substitute(pdIdent(~tmp.zz[, (tmp.num[i] + 
                    1):(tmp.num[i + 1])] - 1), list(i = i)))
                  tmp.block <- c(tmp.block, tmp)
                }
                tmp.block <- lapply(tmp.block, as.formula)
                if (!is.null(grpOrder2)) 
                  tmp.num <- tmp.num[grpOrder2]
                if (!missing(data)) {
                  lmeData <- data
                  lmeData$S <- S
                  lmeData$tmp.zz <- tmp.zz
                  lmeData$tmp.num <- c(tmp.num, rep(0, n - length(tmp.num)))
		  lmeData$tmpGrp<- tmpGrp
                  lme.obj <- lme(formula, random = list(tmpGrp=pdBlocked(tmp.block)), 
                    correlation = correlation, weights = weights, data = lmeData)
                }
                else {
                  lme.obj <- lme(formula, random = list(tmpGrp=pdBlocked(tmp.block)), 
                    correlation = correlation, weights = weights)
                }
            }
            lambda <- as.vector(coef(lme.obj$modelStruct$reStruct))
            lambda <- exp(2 * lambda)
            wei <- NULL
            if (!is.null(lme.obj$modelStruct$corStruct)) {
                wei <- corMatrix(lme.obj$modelStruct$corStruct)
                if (!is.matrix(wei)) {
                  if (length(wei) == 1) 
                    wei <- wei[[1]]
                  else {
                    wei <- bdiag(wei)
                    if (!is.null(grpOrder2)) 
                      wei <- wei[grpOrder2, grpOrder2]
                  }
                }
                cor.est <- lme.obj$modelStruct$corStruct
            }
            if (!is.null(lme.obj$modelStruct$varStruct)) {
                Lambda.wei <- 1/varWeights(lme.obj$modelStruct$varStruct)
                if (!is.null(grpOrder2)) 
                  Lambda.wei <- Lambda.wei[grpOrder2]
                var.est <- lme.obj$modelStruct$varStruct
                if (!is.null(wei)) {
                  wei <- diag(Lambda.wei) %*% wei %*% diag(Lambda.wei)
                }
                else wei <- diag(Lambda.wei^2)
            }
            wei <- (wei + t(wei))/2
            weights <- solve(t(chol(wei)))
            is.lme <- 1
            tmp.Q <- 0
            if (length(lambda) > 1) {
## check the following line
##                for (i in 1:length(Q)) tmp.Q <- tmp.Q + Q[[i]] * lambda[i]
               for (i in 1:length(Q)) tmp.Q <- tmp.Q + Q[[i]] * lambda[i]
                rk.obj <- dsidr(y = y, q = tmp.Q, s = S, weight = weights, 
                  vmu = "m", limnla = 0, tol = ctrl.vals$tol)
               lambda<- 1/lambda
            }
            else {
                rk.obj <- dsidr(y = y, q = Q[[1]], s = S, weight = weights, 
                  vmu = "m", limnla = -log10(lambda[1]), tol = ctrl.vals$tol)
                lambda <- 1/lambda
            }
        }
    }
    else {
        if (length(Q) > 1) {
            rk.obj <- gdmudr(y = y, q = Q, s = S, family = family, 
                varht = varht, vmu = spar, 
                tol1 = ctrl.vals$tol, tol2 = ctrl.vals$tol.g, 
                init = ctrl.vals$init, prec1 = ctrl.vals$prec, 
                prec2 = ctrl.vals$prec.g, maxit1 = ctrl.vals$maxit, 
                maxit2 = ctrl.vals$maxit.g, theta = ctrl.vals$theta)
##            lambda <- 10^rk.obj$theta
	    lambda<- 10^(rk.obj$nlaht-rk.obj$theta)
        }
        else {
            rk.obj <- gdsidr(y = y, q = Q[[1]], s = S, family = family, 
                vmu = spar, maxit = ctrl.vals$maxit.g, varht = varht, 
                limnla = limnla, job = ctrl.vals$job, tol1 = ctrl.vals$tol, 
                tol2 = ctrl.vals$tol.g, prec = ctrl.vals$prec.g)
            lambda <- 10^rk.obj$nlaht
        }
    }
    result <- list(call = match.call(), expand.call = Call, data = dataOrg, 
        y = y, weight = weights, fit = as.vector(rk.obj$fit), 
        s = S, q = Q, lambda = lambda/n, residuals = as.vector(rk.obj$resi), 
        sigma=sqrt(rk.obj$varht), family = family, coef = list(c = rk.obj$c, d = rk.obj$d), 
        df = rk.obj$df, rkpk.obj = rk.obj, scale = scale, cor.est = cor.est, 
        var.est = var.est, control = ctrl.vals)
    if (is.null(S)) {
        result$rkpk.obj$d <- NULL
        result$coef <- list(c = rk.obj$c, d = NULL)
    }
    else {
        result$coef <- list(c = rk.obj$c, d = rk.obj$d)
    }
    class(result) <- "ssr"
    result$call[[1]]<-as.name("ssr")
    if (family != "gaussian") 
        result$weight <- diag(sqrt(result$rkpk.obj$w))
    result$lme.obj <- lme.obj
    result
}
sumList<- 
function(x) 
{ 
## perform summation of components of a list ## 
	if(!is.list(x)) stop("Input must be a list!") 
	val <- 0 
	for(i in 1:length(x)) 
		val <- val + x[[i]] 
	val 
} 

#table2<- 
#function(..., exclude = c(NA, NaN)) 
#{ 
#	args <- list(...) 
#	if((length(args) == 1) && (is.list(args[[1]]))) 
#		args <- args[[1]] 
#	n <- length(args) 
#	if(n == 0) 
#		stop("No arguments") 
#	bin <- 0 
#	lens <- length(args[[1]]) 
#	dims <- numeric(n) 
#	pd <- 1 
#	dn <- vector("list", n) 
#	for(iarg in 1:n) { 
#		i <- args[[iarg]] 
#		if(length(i) != lens) 
#			stop("All arguments must have the same length") 
#		if(!is.category(i)) 
#			i <- category(unlist(i), exclude = exclude) 
#		l <- levels(i) 
#		dims[iarg] <- lenl <- length(l) 
#		dn[[iarg]] <- l 
#		bin <- bin + pd * (as.numeric(i) - 1) 
#		pd <- pd * lenl 
#	} 
#	names(dn) <- names(args) 
#	if(length(wna <- which.na(bin))) 
#		bin <- bin[ - wna] 
#	array(tabulate(bin + 1, pd), dims, dn) 
#} 
tp<-
function(s, u = s, order = 2)
{
## calculate true thin-plate spline kernel
	tp.p <- tp.pseudo(s = s, u = u, order = order)
	swk <- tp.term(s, order = order)
	qrmat <- qr(swk)
	smat <- qr.Q(qrmat, complete = TRUE)
	qmat <- smat[, (ncol(swk) + 1):(nrow(swk))]
	smat <- smat[, 1:(ncol(swk))]
	if(!missing(u)) {
		twk <- tp.term(u, order = order)
		qrmat2 <- qr(twk)
		smat2 <- qr.Q(qrmat2, complete = TRUE)
		qmat2 <- smat2[, (ncol(twk) + 1):(nrow(twk))]
		smat2 <- smat2[, 1:(ncol(twk))]
	}
	else {
		qmat2 <- qmat
		smat2 <- smat
	}
	qmat <- qmat %*% t(qmat) %*% tp.p %*% qmat2 %*% t(qmat2)
	qmat <- diag(smat[, 1]) %*% qmat %*% diag(smat2[, 1])	
	qmat
}
tp.pseudo<-
function(s, u = s, order = 2)
{
## Calculate pseudo thin-plate spline kernel
        
	if(is.matrix(s)) {
		d <- ncol(s)
		s <- unlist(as.vector(s))
	}
	else {
		if(is.list(s)) {
			d <- length(s)
			s <- unlist(s)
		}
                else d<- 1
	}
	if(missing(u))
		u <- s
	else {
		if(is.matrix(u)) {
			if(ncol(u) != d)
				stop("u must match")
			else u <- unlist(as.vector(u))
		}
		else {
			if(is.list(u)) {
				if(length(u) != d)
				  stop("u must match")
				else u <- unlist(u)
			}
		}
	}
	r <- 2 * order - d
	if(r <= 0)
		stop(" positive 2*order-d is required")
	N <- length(s)/d
	M <- length(u)/d
	resul <- matrix(.C("tp_ker",
		as.double(s),
		as.double(u),
		as.integer(d),
		as.integer(order),
		as.integer(N),
		as.integer(M),
		val = double(N * M),
		PACKAGE="assist")$val, ncol = M, byrow = TRUE)
	if((d %% 2) == 0)
		resul * ((-1)^(d/2 + order + 1))
	else resul * sign(sin(pi * (d/2 - order)))
}
tp.term<-
function(x, order = 2)
{
## calculate bases for thin-plate splines
	if(!is.list(x))
		x <- as.matrix(x)
	if(!is.matrix(x)) {
		tmp <- NULL
		for(i in 1:length(x)) {
			tmp <- cbind(tmp, as.vector(x[[i]]))
		}
		x <- tmp
	}
	d <- ncol(x)
	result <- .C("tp_term",
		as.integer(d),
		as.integer(order),
		p = integer(order^d * d),
		PACKAGE="assist")$p
	result <- matrix(result, nrow = d, byrow = FALSE)
	term <- choose(order + d - 1, d)
	mat <- NULL
	for(i in 1:term)
		mat <- rbind(mat, apply(t(x)^result[, i], 2, prod))
	t(mat)
}
tp.linear<-function(s, u=s){
    if(is.vector(s)) return(kron(s,u))
    if(is.list(s)){
		d<- length(u)
		s<-matrix(unlist(s), ncol=d, byrow=F) 
	}
	d<- ncol(s)
     Tx <- cbind(1,s)
     Rx <- qr.R(qr(Tx))
     zx <- Tx%*%solve(Rx)
     xx <- zx[,-1]	 
	 if(missing (u)) sumList(sapply(1:d, function(a, b) list(kron(b[,a])), b=xx))
	 else{
	 	if(is.list(u)){
			dy<- length(u)
			if(dy!=d) stop("u must have the same length as s")
		    u<-matrix(unlist(u), ncol=dy, byrow=F) 
	    }
	    dy<- ncol(u)
     	Ty <- cbind(1,u)
     	Ry <- qr.R(qr(Ty))
     	zy <- Ty%*%solve(Ry)
     	yy <- zy[,-1]
     	sumList(sapply(1:d, function(a, b, f) list(kron(b[,a], f[,a])), b=xx, f=yy))
    }
}

xyplot2<-
function(formula, data, type = "l", ...)
{
## this is to extend xyplot, and multiple datasets input ##
## are allowed                                           ##
## further work is needed for completion                 ##
	call <- match.call()
#	var.names <- c(as.character(call$formula[[2]]), as.character(call$
#		formula[[3]][[2]]), as.character(call$formula[[3]][[3]]))
        var.names<- all.vars(call$formula)
	if(is.data.frame(data)) {
		plot.resu <- xyplot(formula, data = data)
	}
	else {
		tmp.dat <- data[[1]][var.names]
		index.names <- rep("dat1", nrow(data[[1]]))
		for(i in 2:length(data)) {
			tmp.dat <- rbind(tmp.dat, data[[i]][var.names])
			index.names <- c(index.names, rep(paste("dat", i, sep
				 = ""), nrow(data[[i]])))
		}
		tmp.dat$index.names <- index.names
		if(missing(type))
			type <- rep("l", length(data))
		plot.resu <- xyplot(formula, data = tmp.dat, subscripts = TRUE, 
			groups = index.names, strip = function(...)
		strip.default(..., style = 1), ..., panel = function(x, y, 
			subscripts, groups, ...)
		{
			dat1 <- groups[subscripts] == "dat1"
			panel.xyplot(x[dat1], y[dat1], type = type[1])
			ind <- unique(groups)
			for(i in 2:length(ind)) {
				fit <- groups[subscripts] == ind[i]
				panel.superpose(x[fit], y[fit], subscripts[fit],
				  groups, type = type[i], lty = 3)
			}
		}
		)
	}
	plot.resu
}
prodDiag<-
function(q, x)
{
	t(t(q) * x)
}

print.ssr<- 
function(x, ...) 
{ 
	cat("Smoothing spline regression fit by ") 
	switch(x$rkpk.obj$vmu, 
		v = cat("GCV"), 
		m = cat("GML"), 
		u = cat("UBR"), 
		"~u" = cat("~UBR")) 
	cat(" method\n") 
	if(!is.null(cl <- x$call)) { 
		cat("Call: ") 
		dput(cl) 
	} 
	switch(x$rkpk.obj$vmu, 
		v = cat("\nGCV"), 
		m = cat("\nGML"), 
		u = cat("\nUBR"), 
		"~u" = cat("\n~UBR")) 
#	if(length(x$lambda) == 1) 
#		cat(" estimate(s) of smoothing parameter(s) :", format(x$lambda 
#			), "\n") 
#	else cat(" estimate(s) of smoothing parameter(s) :", format(1/x$lambda), 
#			"\n") 
	cat(" estimate(s) of smoothing parameter(s) :", format(x$lambda 
			), "\n")
	cat("Equivalent Degrees of Freedom (DF): ", format(x$rkpk.obj$df), "\n" 
		) 
	cat("Estimate of sigma: ", format(sqrt(x$rkpk.obj$varht)), "\n") 
	if(!is.null(x$cor.est)) 
		print(x$cor.est) 
	if(!is.null(x$var.est)) 
		print(x$var.est) 
	cat("\nNumber of Observations: ") 
	if(is.matrix(x$y)) 
		cat(length(x$y[, 1]), "\n") 
	else cat(length(x$y), "\n") 
} 

predict.ssr<- function(object, newdata=NULL,terms, pstd=TRUE, ...)
{
## calculate posterior STD and fits 
## of SS ANOVA components from ssr object
    Call<- match.call()
    if(length(object$q)>1) is.dmudr<-TRUE
     else is.dmudr<-FALSE

    if(inherits(object$cor.est, "corStruct") || inherits(object$var.est, "varFunc"))
         lmeCalled<- TRUE
    else lmeCalled<- FALSE
    nobs<- object$rkpk.obj$nobs
    nnull<- object$rkpk.obj$nnull
    if(length(object$lambda)==1) theta<-1
#    else theta<- object$lambda*nobs
    else theta<- 1/(object$lambda*nobs)

    if(is.null(newdata)) eval.len<- nobs
    else   eval.len<- length(newdata[[1]])

## if scale=TRUE
    if(object$scale==TRUE && !is.null(newdata)){
        varName<- intersect(unique(all.vars(object$call$rk)), names(newdata))
        if(!is.data.frame(object$data)) {        
         data<- as.list(varName)
         data<- data.frame(lapply(data, get))
         names(data)<-varName
              }	
         else data<- object$data
         for(var in varName){
            newdata[[var]]<- ident(data[[var]], newdata[[var]])
            data[[var]]<- ident(data[[var]])              
            }
       data.used<- data
   }
   else data.used<- object$data

## begin-of-old-scaling 
## if scale=TRUE 
#     if(object$scale==TRUE){
#        if(!is.null(newdata)) newdata<- ident(data.frame(object$data), newdata)
#        data.used<- data.frame(ident(object$data))
#      }
#      else  data.used<- object$data
## end-of-old-scaling

    if(missing(terms) || is.null(terms)) terms<- rep(1, nnull+ length(theta))
    terms<-as.matrix(terms)

    if((ncol(terms) != nnull+length(theta)) && (nrow(terms) != nnull+length(theta)))
       stop(" the input of terms must match ") 
   
    if(ncol(terms)!= nnull+length(theta)) terms<-t(terms)

   if (nnull > 0) {
        terms1 <- matrix(terms[, 1:(ncol(terms) - length(theta))], 
            nrow = nrow(terms))
        terms2 <- matrix(terms[, (ncol(terms1) + 1):(ncol(terms))], 
            nrow = nrow(terms))
    }
    else {
        terms1 <- NULL
        terms2 <- terms
    }


# calculate S and Q and r's
    swk<- object$s
    qwk<- object$q
    if(is.null(newdata)){
           phi<- object$s
           rwk<- Rwk<- object$q
          } 
    else{ 
    rwk<- rkEval(object$expand.call$rk, data.used, newdata)
## if subset
        if(!is.null(eval(object$call$subset, envir=object$data))){
           for(compon in 1:length(rwk)){
              rwk[[compon]]<- rwk[[compon]][eval(object$call$subset, envir=object$data),]
            }
         }

         if(nnull>0){
		phi<- model.matrix(eval(object$expand.call$formula[-2]), data=newdata)

        	if(nrow(phi)==1) 
             	phi<- matrix(rep(as.vector(phi),dim(newdata)[1]), nrow=dim(newdata)[1])
        } 
        Rwk<- eval(object$call$rk, newdata)
        if(is.matrix(Rwk)) Rwk<- list(Rwk)
     }
  
    phi.new<- NULL
    if(nnull>0){
    	for(i in 1:ncol(swk)){
      		phi.new<- cbind(phi.new, as.vector(phi[,i] %*% t(as.matrix(terms1[, i]))))
     	}
    }
   if(is.matrix(rwk)) rwk<- list(rwk) 
   if(is.matrix(Rwk)) Rwk<- list(Rwk)                  

   qSum<- Rsum<- 0
   xi<- xi2<- 0
    for(i in 1:length(theta)){
      rwk[[i]]<-  rwk[[i]]*theta[i]
      xi2<- xi2+ as.vector(rwk[[i]]) %*% t(as.matrix(terms2[, i]))
      Rsum<- Rsum+ (as.vector(Rwk[[i]])*theta[i]) %*% t(as.matrix(terms2[, i]))
      qSum<- qSum+ qwk[[i]]*theta[i]

      if(!is.null(object$weight)) rwk[[i]]<-object$weight%*%rwk[[i]]
      xi<- xi+ as.vector(rwk[[i]]) %*% t(as.matrix(terms2[, i]))
     }


# call dsidr
   y<-object$y
   if(is.dmudr && (lmeCalled==FALSE)){
      if(is.null(object$expand.call$family) || 
              object$expand.call$family == "gaussian"){
           dsidr.fit<- dsidr(s=swk, q=qSum, y=y, 
                    weight=object$weight,
                    vmu=object$call$spar, 
	            varht=object$expand.call$varht, 
                    tol=object$control$tol,
                    job=object$control$job, 
                    limnla=object$call$limnla)
              }
      else{
           dsidr.fit<- gdsidr(s=swk, q=qSum, y=y, 
              family=object$expand.call$family,
              vmu=object$expand.call$spar,
	      varht=object$expand.call$varht, 
              tol1=object$control$tol,
              tol2=object$control$tol.g, 
              maxit=object$control$maxit.g,
              job=object$control$job, 
              limnla=object$expand.call$limnla,
              prec=object$control$prec)
     }
   }
   else dsidr.fit<- object$rkpk.obj

## pstd=TRUE   
   if(pstd){
      b<- dsidr.fit$varht/(10^dsidr.fit$nlaht)

# call dsms and dcrdr   
    if(nnull>0) dsmsV<- matrix(dsms(dsidr.fit)$sms, ncol=nnull)
     cr<-0
     dr<-0

     for(i in 1:length(rwk)){
       dcrdrV<- dcrdr(dsidr.fit, rwk[[i]])
       cr<- cr + dcrdrV$cr %*% t(as.matrix(terms2[, i]))
       dr<- dr + dcrdrV$dr %*% t(as.matrix(terms2[, i]))
       }
# collect STDs
#  variance for the rk part
      ciRK<-  matrix( as.vector(cr) *as.vector(xi) , ncol=nobs, byrow=TRUE)
      ciRK<- matrix(apply(ciRK, 1, sum), nrow=eval.len, byrow=FALSE)
      Rsum<- apply(Rsum, 2, function(x, r) x[seq(1, length(x), length=r)], 
             eval.len)
      ciRK<- Rsum-ciRK

      if(nnull>0){
#  variance for the  main effect part
      	ciMain<- apply(phi, 1, function(x, terms1, w) diag(terms1 %*%(x*w)%*% (x*t(terms1))), terms1=terms1, w=dsmsV)
        ciMain<- matrix(as.vector(ciMain), nrow=nrow(terms), byrow=FALSE)

# covariance between both parts
        ciCross<-  matrix(as.vector(t(phi.new)) * as.vector(dr), ncol=nnull, byrow=TRUE)
        ciCross<-  matrix(apply(ciCross, 1, sum), nrow=eval.len, byrow=FALSE)
# overall covariance
        ci<- t(ciMain) + ciRK -2* ciCross
      }
      else ci<- ciRK
     
      ci<- ci * (ci>0)
    }
# caculate fits
     ccc<- dsidr.fit$c
     if(nnull>0){
     	fit<- phi %*% matDiag.prod(t(terms1), dsidr.fit$d) + apply(xi2, 2, 
       		function(x, w, r) matVec.prod(matrix(x, ncol=r, byrow=TRUE), as.vector(w), left=FALSE),w=ccc, r=nobs)
     }
     else{
        fit<-  apply(xi2, 2, function(x, w, r) matVec.prod(matrix(x, ncol=r, byrow=TRUE), as.vector(w), left=FALSE),w=ccc, r=nobs)
     }
          
     if(ncol(fit)==1) fit<- as.vector(fit)
     if(pstd){
        if(ncol(ci)==1) ci<- as.vector(ci)
        resul<-list(fit=fit, pstd=sqrt(ci*b))
       }
     else resul<-list(fit=fit, pstd=NULL)
   class(resul)<-c("bCI", "list")
   resul
}
    
summary.ssr<- function(object, ...){
##summary facility for ssr object
   resul<-list(call=object$call,
                sigma=sqrt(object$rkpk.obj$varht), 
                spar=object$rkpk.obj$vmu,
                lambda=object$lambda, 
                df=object$df,
                nobs=nrow(object$s),
                cor.est=object$cor.est,
                var.est=object$var.est,
                coef.c=object$coef$c,               
                coef.d=object$coef$d)
    resul$family<- "gaussian"
    if(!is.null(d.name<-dimnames(object$s)[[2]])) names(resul$coef.d)<- d.name   
    if(!is.null(object$call$family)) resul$family<- object$call$family
    class(resul)<- "summary.ssr"
    resul
}
print.summary.ssr<- function(x, ...){  
       cat("Smoothing spline regression fit by ")
       switch(x$spar,
		"v"= cat("GCV"),
		"m"= cat("GML"),
		"u"= cat("UBR"),
		"~u"= cat("~UBR"))
       cat(" method\n")

	if(!is.null(cl <- x$call)) {
		cat("Call: ")
		dput(cl)
	}
        cat("\nCoefficients (d):\n")
        print(x$coef.d)

      switch(x$spar,
		"v"= cat("\nGCV"),
		"m"= cat("\nGML"),
		"u"= cat("\nUBR"),
		"~u"= cat("\n~UBR"))
#	if(length(x$lambda)==1)
#               cat(" estimate(s) of smoothing parameter(s) :", format(x$lambda), "\n")
#        else   cat(" estimate(s) of smoothing parameter(s) :", format(1.0/x$lambda), "\n")
	cat(" estimate(s) of smoothing parameter(s) :", format(x$lambda), "\n")
	cat("Equivalent Degrees of Freedom (DF): ", format(x$df), "\n")
	cat("Estimate of sigma: ", format(x$sigma), "\n")
       if(!is.null(x$cor.est)) print(x$cor.est)
        if(!is.null(x$var.est)) print(x$var.est)
       cat("\nNumber of Observations: ", x$nobs, "\n")
  }   

print.anova.ssr<-function(x, ...){
    cat("\nTesting H_0: f in the NULL space\n")
    cat("\n")

    if(!is.null(x$lmp.test)){
       LMP<- c(format(x$lmp.test[[1]]), x$simu.size, 
             format(mean(x$lmp.test[[2]]> x$lmp.test[[1]])))    
      }  
    if(!is.null(x$gcv.test)){       
       GCV<- c(format(x$gcv.test[[1]]), x$simu.size,
               format(mean(x$gcv.test[[2]]< x$gcv.test[[1]])))      
       x<-data.frame(rbind(LMP, GCV))       
       names(x)<-c("test.value", "simu.size", "simu.p-value")
     }
   if(!is.null(x$gml.test)){       
       pp<- 0.5-pchisq(-2*log(x$gml.test[[1]])*(x$dim[1]-x$dim[2]), 1)*0.5            
       GML<-c(format(x$gml.test[[1]]), x$simu.size, 
                format(mean(x$gml.test[[2]]< x$gml.test[[1]])))
       x<- data.frame(rbind(LMP, GML))
       x<- cbind(x, c("",as.character(format(pp)))) 
       row.names(x)<- c("LMP", "GML")
       names(x)<-c("test.value", "simu.size", "simu.p-value", "approximate.p-value")
     }
   print(x)
}

anova.ssr<- function(object, simu.size=100, ...){
## anova only for gaussion single smoothing parameters
## calculate testing statistics
 if((object$family=="gaussian") && (length(object$q)==1) && is.null(object$cor.est) && is.null(object$var.est)){ 
    nobs<- length(object$y)
    nobs.par<- ncol(object$s)
    qr.result<- qr(object$s)
    qrq<- qr.Q(qr.result, complete=TRUE)
    qrq2<- qrq[, (ncol(object$s)+1): nobs]
    yy<- t(qrq2) %*% object$y
    u<- t(qrq2)%*%object$q[[1]] %*% qrq2
    l<- eigen((u+t(u))/2)
    yy<- as.vector(t(l$vectors) %*% yy)
    lambda<- l$values
    rv<- 1+ lambda/(10^object$rkpk.obj$nlaht)
    tt1<- sum(lambda* (yy^2))/sum(yy^2)
## using simulation to calculate p-values
    z<-matrix(rnorm(simu.size*length(rv)), ncol=simu.size)*sqrt(object$rkpk.obj$varht)
    if(object$rkpk.obj$vmu=="v"){
            tt2<- sum((yy/rv)^2)/sum(1/(rv^2))/sum(yy^2)
            z2<-apply(z, 2, 
                  function(x, a) sum((x/a)^2)/sum(1/(a^2))/sum(x^2), a=rv)
            gcv.test<-list(value=tt2, simu=z2)
            gml.test<- NULL
            }
    else{
        if(object$rkpk.obj$vmu=="m"){
            tt2<- sum(yy^2/rv)*exp(mean(log(rv)))/sum(yy^2)
            z2<- apply(z, 2,
                  function(x, a) sum(x^2/a)*exp(mean(log(a)))/sum(x^2), a=rv)
            gml.test<-list(value=tt2, simu=z2)
            gcv.test<- NULL
          }
        }
    z1<- apply(z, 2, function(x, a) sum(a*x^2)/sum(x^2), a=lambda)
    lmp.test<-list(value=tt1, simu=z1)
    resul<-list(lmp.test=lmp.test, gml.test=gml.test,  gcv.test=gcv.test, 
              dim=c(nobs, nobs.par), simu.size=simu.size)
    class(resul)<-"anova.ssr"
    resul
   }
 else stop("ANOVA can not be calculated!")
}


deviance.ssr <-
function (object, residuals = FALSE, ...) 
{
    if (is.null(object$family)) 
        object$family <- "gaussian"
    switch(f <- object$family, gaussian = {
        resu <- object$y - object$fit
        if (!is.null(object$weight)) 
            resu <- as.vector(object$weight %*% resu)
        if (!residuals) 
            sum(resu^2)
        else resu
    }, binary = {
        resu<-binomial()$dev.resids(mu=object$fit, y=object$y, wt=diag(object$weight))
	if(residuals) resu
	else sum(resu^2)
    }, binomial = {
        resu<- binomial()$dev.resids(mu=object$fit, y=object$y, wt=diag(object$weight))
	if(residuals) resu
	else sum(resu^2)
    }, poisson = {
        resu<-poisson()$dev.resids(mu=object$fit, y=object$y, wt=diag(object$weight))
	if(residuals) resu
	else sum(resu^2)
    }, gamma = {
        resu<-Gamma()$dev.resids(mu=object$fit, y=object$y, wt=diag(object$weight))
	if(residuals) resu
	else sum(resu^2)
    })
}

# removed since it creates a NOTE
#.First.lib <- function(lib, pkg)
#{
#    library.dynam("assist", pkg, lib)
#}
#intervals<- function(object, ...) UseMethod("intervals")

#### NNR-Part
nnr<- function(formula, func, spar="v", 
               data=sys.parent(), start=list(), verbose=FALSE,  control=list()) 
{ 
## fit a general nonlinear spline model 
## using backfitting algorithm 
  thisCall <- match.call() 
## figure out the response 
  y<- eval(getResponseFormula(formula)[[2]], data) 
  nobs<- length(y) 
 
## get info for "f" 
  f.info<-getFunInfo(func) 
  funcName<- f.info$fName 
  funcArg<- f.info$f.arg 
  thisCall$rk<- f.info$fRkForm 
  thisCall$formF<- f.info$fNullForm 
  lengthF<- length(funcName) 
  thisCall$formF[[lengthF+1]]<- 1
## construct $f$ 
  splineF<- splineF0<- list(lengthF) 
  for(numF in 1:lengthF){ 
     splineF[[numF]]<-function(...){ 
       thisFunName<- as.character(match.call()[[1]]) 
       thisFunOrder<-match(thisFunName, funcName) 
       fEst[[thisFunOrder]] 
     } 
             
     splineF0[[numF]]<-function(...) { 
       thisFunName<- as.character(match.call()[[1]]) 
        thisFunOrder<-match(thisFunName, funcName) 
       res<-data.frame(cbind(...)) 
       tmp.delta3<- delta3 
       tmp.delta3[[thisFunOrder]]<- res 
       assign("delta3", tmp.delta3, envir=environment(eval(parse(text=thisFunName)))) 
#       res 
       1
     } 
   }    
 
## get args of f 
   for(i in 1:lengthF){ 
     assign(funcName[i], splineF0[[i]]) 
    } 
   delta3<- list() 
   if(missing(data)) gb<-eval(formula[[3]]) 
   else gb<-eval(formula[[3]], data) 
   Smat<-Qmat<- list(lengthF) 
   Smat[[lengthF+1]]<- 0
   for(num in 1:lengthF){ 
      names(delta3[[num]])<-funcArg[[num]] 
      if(length(thisCall$formF)>0 && !is.null(thisCall$formF[[num]])) 
           Smat[[num]]<- model.matrix2(thisCall$formF[[num]], delta3[[num]]) 
      Qmat[[num]]<- eval(thisCall$rk[[num]], delta3[[num]]) 
   } 
 
## get start values for f 
   if(!missing(data)) tmpData<- data.frame(delta3, data) 
   else tmpData<- delta3 
   fhat0<- eval(thisCall$start, tmpData) 
   if(length(fhat0)!=lengthF) stop("Invalid input of start values") 
   if(!is.null(names(fhat0))) fhat0<- fhat0[funcName]
   fhat<-fhat0<- lapply(fhat0, function(x,y) if(length(x)>1) x else rep(x,y), y=nobs) 
 
## get control info 
   defCtrl<- nnr.control()  
   if(!missing(control) && !is.null(control)) 
     for(nam in names(control)) defCtrl[[nam]]<-control[[nam]] 
   
   defCtrl$vmu<-spar 
    
## iteration starts  
   iteration<-0 
   rss<- 1 
   rssOR<- c() 
   incDelta<- defCtrl$increment 
   rkFit<- wList<- yList<- list() 
   pls<- NULL 
   repInd<- rep(1:lengthF, rep(defCtrl$backfit, lengthF)) 
 
   repeat{       
      iteration<- iteration+1       
      if(verbose) cat("\nIteration ", iteration, "\n") 
      for(i in 1:lengthF){ 
           assign(funcName[i], splineF[[i]]) 
        } 
 
## begin backfitting 
      for(numF in repInd){  
        fEst<- lapply(fhat, function(x) cbind(x, x, x))      
## calculate D and E  
        fEst[[numF]]<- cbind(fhat[[numF]], fhat[[numF]]+incDelta, fhat[[numF]]-incDelta) 
 
        if (missing(data)) fstD<- eval(formula[[3]]) 
        else  fstD<- eval(formula[[3]], data) 
        Dmat<- (fstD[,2]-fstD[,3])/incDelta/2.0 
         
        if(defCtrl$converg=="ortho"){ 
            rssOR[numF]<- abs(sum((y-fstD[,1])*Dmat))/sqrt(sum((y-fstD[,1])^2)*sum(Dmat^2)) 
        } 
  
        if(defCtrl$method=="GN") Emat<-0 
        else{ 
           sndD<- (fstD[,2]+fstD[,3]-2*fstD[,1])/incDelta/incDelta 
           Emat<- sndD*(y-fstD[,1]) 
          } 
 
        DeltaMat<- abs(Dmat^2-Emat) 
        uVec<- Dmat*(fstD[,1]-y) 
        yTilde<-fhat[[numF]]-uVec/DeltaMat 
 
## call dsidr/dmudr 
        if(!is.list(Qmat[[numF]])) Qmat[[numF]]<- list(Qmat[[numF]]) 
        if(length(Qmat[[numF]])==1) { 
            rkFit[[numF]]<- dsidr(y=yTilde, q=Qmat[[numF]][[1]],  
                   s=Smat[[numF]], weight=diag(sqrt(DeltaMat)), 
                   tol=defCtrl$tol, vmu=defCtrl$vmu, 
                   job=defCtrl$job,  limnla=defCtrl$limnla, 
		   varht=defCtrl$varht) 
 
         } 
         else{ 
            rkFit[[numF]]<- dmudr(y=yTilde, q=Qmat[[numF]],  
                   s=Smat[[numF]], weight=diag(sqrt(DeltaMat)), 
                   tol=defCtrl$tol,  vmu=defCtrl$vmu, 
                   prec=defCtrl$prec, maxit=defCtrl$maxit, 
                   init=defCtrl$init, varht=defCtrl$varht, 
		   theta=defCtrl$theta) 
         } 
 
    if(verbose) cat("Estimate d: ", format(rkFit[[numF]]$d), "\n") 
     fhat[[numF]]<- rkFit[[numF]]$fit 
     wList[[numF]]<- diag(sqrt(DeltaMat)) 
     yList[[numF]]<- yTilde 
     } 
    if(defCtrl$converg =="ortho") rss<- max(rssOR) 
    else{ 
        rss<- sum(abs(unlist(fhat)-unlist(fhat0)))/(1+sum(abs(unlist(fhat0))))               
    } 
      
    if(verbose) cat("Convergence Criterion: ", format(rss), "\n") 
    if((rss<defCtrl$toler) || (iteration>=defCtrl$max.iter)) break  
    fhat0<- fhat      
  } 
 
   if(iteration>defCtrl$max.iter) warnings("Convergence not reached!") 
   forCI<- list(expand.call=thisCall, rkpk.obj=rkFit, y=yList,  
         s=Smat, q=Qmat, data=delta3, control=defCtrl, weight=wList) 
   df<- unlist(lapply(rkFit, function(x) x$df)) 
#   lambda<- unlist(lapply(rkFit, function(x) ifelse(is.null(x$theta),x$nlaht,-1*x$theta)))
   lambda<- unlist(lapply(rkFit, function(x) ifelse(is.null(x$theta),x$nlaht,x$nlaht-x$theta))) 
   names(fhat)<- funcName 
   result<-list(call=match.call(),  
        data=data, 
        funcFitted=fhat,  
        df=list(total=nobs, f=sum(df)), 
        lambda=lambda, 
        forCI=forCI,            
        iteration=list(times=iteration, rss=rss, maxIter=defCtrl$max.iter)) 
  result$residuals<- y-fstD[,1] 
  result$sigma<-sqrt(sum(result$residuals^2)/(nobs-result$df$f)) 
  class(result)<- c("nnr") 
  result 
} 
       
nnr.control<- function(job = -1, tol = 0, max.iter=50, 
             init = 0, limnla=c(-10, 0), varht=NULL, theta= NULL,  
             prec = 1e-006, maxit = 30, method="NR", increment=0.0001,  
             backfit=5, converg="coef", toler=0.001) 
{ 
  if((init == 1) && is.null(theta)) 
	stop("initial values for theta is needed") 
  if(!(converg =="ortho" || converg=="coef")) stop("unknown convergence approach!") 
  if(!(method=="NR" || method=="GN")) stop("unknown approximation method!") 
  list(job = job, tol = tol, max.iter=max.iter,init = init, limnla=limnla,  
      varht=varht, theta= theta, prec = prec, maxit = maxit, method=method, increment=increment,  
       backfit=backfit, converg=converg, toler=toler) 
} 

print.nnr<-
function(x, ...)
{
   if(inherits(x, "nnr")){
     cat("Nonlinear Nonparametric Regression Model Fit by ")
     if(x$forCI$control$method=="GN") cat("Gauss-Newton Method\n")
     else cat("Newton-Raphson Method\n")
     }

   cat(" Model:", deparse(as.vector(x$call$formula)), "\n")
   if(!is.null(x$call$data))
     cat(" Data:", deparse(x$call$data), "\n")
   cat("\n")   
     switch(x$forCI$rkpk.obj[[1]]$vmu,
          "v"= cat(" GCV"),
          "m"= cat(" GML"),
          "u"= cat(" UBR"))
  cat(" estimate(s) of smoothing parameter(s):", format(10^x$lambda/x$df$total), "\n")
  cat(" Equivalent Degrees of Freedom (DF): ", format(sum(x$df$f)), "\n")
  cat("\n Residual standard error:", format(x$sigma), "\n")
  cat(" Number of Observations:", length(x$forCI$y[[1]]), "\n") 
  if(x$iter$times<x$iter$maxIter)
        cat(" Converged after", format(x$iter[1]), "iterations\n")
  else cat(" Not Converged after", format(x$iter[1]), "iterations\n")
}

print.summary.nnr<- function(x, ...){     
   cat("Nonlinear Nonparametric Regression Model Fit by ")
   if(x$GN) cat("Gauss-Newton Method\n")
   else cat("Newton-Raphson Method\n")     
   cat("  Model:", deparse(as.vector(x$call$formula)), "\n")
   if(!is.null(x$call$data)) cat("  Data:", deparse(x$call$data), "\n")
   cat("\n")

   switch(x$vmu,
          "v"= cat("GCV"),
          "m"= cat("GML"),
          "u"= cat("UBR"))
  cat(" estimate(s) of smoothing spline parameter(s):", format(x$lamb), "\n")  
  cat("Equivalent Degrees of Freedom (DF) for spline function: ", format(x$df), "\n")
  cat("Residual standard error:", format(x$sigma), "\n")
  if(x$iter$times<x$iter$maxIter)
        cat("\nConverged after", format(x$iter[1]), "iterations\n")
  else cat("\nNot Converged after", format(x$iter[1]), "iterations\n")
}

summary.nnr<- function(object, ...){
## summary for nnr fit   ##
  resul<-list(call=object$call, 
          sigma=object$sigma,
          iter=object$iteration)
  resul$vmu<- object$forCI$rkpk.obj[[1]]$vmu
  spar<- NULL
  for(i in 1:length(object$forCI$rkpk.obj)){
     if(object$forCI$rkpk.obj[[i]]$nq>1) spar<- c(spar, 10^(object$forCI$rkpk.obj[[i]]$nlaht-object$forCI$rkpk.obj[[i]]$theta))
     else spar<- c(spar, 10^object$forCI$rkpk.obj[[i]]$nlaht)
     }
  resul$lamb<- spar/object$df$total
  resul$df<- object$df$f
  resul$GN<-object$forCI$control$method=="GN"      
  class(resul)<- "summary.nnr"
  resul
}

intervals.nnr<- function(object, level=0.95, newdata=NULL, terms, pstd=TRUE, ...) 
{ 
## prediction, bayesian CI for nnr fit 
   Call<- match.call() 
   if(!inherits(object, "nnr")) stop("Input doesn't match") 
   resul<-list() 
 
## spline structures 
   f.info<-getFunInfo(eval(object$call$func)) 
   funcName<- f.info$fName 
   funcArg<- f.info$f.arg 
   lengthF<- length(funcName) 
 
## get newdata 
   if(!missing(newdata)){ 
     splineF0<- list(lengthF) 
     for(numF in 1:lengthF){             
       splineF0[[numF]]<-function(...) { 
          thisFunName<- as.character(match.call()[[1]]) 
          thisFunOrder<-match(thisFunName, funcName) 
          res<-data.frame(cbind(...)) 
          tmp.delta3<- delta3 
          tmp.delta3[[thisFunOrder]]<- res 
          assign("delta3", tmp.delta3, envir=environment(eval(parse(text=thisFunName)))) 
          res 
          } 
       }     
 
     delta3<- list() 
     for(i in 1:lengthF){ 
     assign(funcName[i], splineF0[[i]]) 
     } 
     eval(object$call$formula[[3]], newdata)   
   } 
   else delta3<- object$forCI$data 
## calculate CI one by one 
   for(numF in 1:lengthF){ 
      ssr.obj<- object$forCI 
      ssr.obj$rkpk.obj<- ssr.obj$rkpk.obj[[numF]] 
      ssr.obj$y<- ssr.obj$y[[numF]] 
      ssr.obj$weight<- ssr.obj$weight[[numF]] 
      ssr.obj$s<- ssr.obj$s[[numF]] 
      ssr.obj$q<- ssr.obj$q[[numF]] 
      newdata<- delta3[[numF]] 
      if(!is.null(newdata))  names(newdata)<- funcArg[[numF]] 
 
      if(length(ssr.obj$q)>1) is.dmudr<-TRUE 
      else is.dmudr<-FALSE 
      if(is.null(theta<- ssr.obj$rkpk.obj$theta)) theta<- 1 
      else theta<- 10^theta  
 
 
      nobs<- length(ssr.obj$y) 
      nnull<- ssr.obj$rkpk.obj$nnull 
 
      if(is.null(newdata)) eval.len<- nobs 
      else   eval.len<- nrow(newdata) 
 
      if(!missing(terms)) cTerm<- terms[[funcName[numF]]] 
      else cTerm<- NULL 
 
      if(is.null(cTerm)) cTerm<-  rep(1, nnull+ length(theta)) 
      cTerm<- as.matrix(cTerm) 
 
      if((ncol(cTerm) != nnull+length(theta)) && (nrow(cTerm) != nnull+length(theta))) 
           stop(" the input of terms must match ")  
    
      if(ncol(cTerm)!= nnull+length(theta)) cTerm<-t(cTerm) 
      if(nnull>0){ 
        terms1<- matrix(cTerm[, 1:(ncol(cTerm)-length(theta))], nrow=nrow(cTerm)) 
        terms2<- matrix(cTerm[, (ncol(terms1)+1):(ncol(cTerm))], nrow=nrow(cTerm)) 
      } 
      else{ 
        terms1<- NULL 
        terms2<- cTerm 
      } 
 
      swk<- ssr.obj$s 
      qwk<- ssr.obj$q 
 
      if(is.null(newdata)){ 
           phi<- ssr.obj$s 
           rwk<- Rwk<- ssr.obj$q 
          } 
      else{        
          phi<- NULL      
          if(!is.null(ssr.obj$expand.call$formF[[numF]])){ 
              phi<- model.matrix2(ssr.obj$expand.call$formF[[numF]], newdata)              
          } 
          Rwk<- eval(ssr.obj$expand.call$rk[[numF]], newdata) 
          if(!is.list(Rwk)) Rwk<- list(Rwk)                   
          rwk<- rkEval(ssr.obj$expand.call$rk[[numF]], ssr.obj$data[[numF]], newdata) 
          if(!is.list(rwk)) rwk<- list(rwk)           
      } 
 
      phi.new<- NULL 
      if(!is.null(terms1)){ 
         for(i in 1:ncol(swk)){ 
            phi.new<- cbind(phi.new, as.vector(phi[,i] %*% t(as.matrix(terms1[, i])))) 
         } 
      } 
      if(is.matrix(rwk)) rwk<- list(rwk)  
      if(is.matrix(Rwk)) Rwk<- list(Rwk)                   
      qSum<- Rsum<- 0 
      xi<- xi2<- 0 
      for(i in 1:length(theta)){ 
         rwk[[i]]<-  rwk[[i]]*theta[i] 
         xi2<- xi2+ as.vector(rwk[[i]]) %*% t(as.matrix(terms2[, i])) 
         Rsum<- Rsum+ (as.vector(Rwk[[i]])*theta[i]) %*% t(as.matrix(terms2[, i])) 
         qSum<- qSum+ qwk[[i]]*theta[i] 
 
         if(!is.null(ssr.obj$weight)) rwk[[i]]<-ssr.obj$weight%*%rwk[[i]] 
         xi<- xi+ as.vector(rwk[[i]]) %*% t(as.matrix(terms2[, i])) 
     } 
# call dsidr 
     y<-ssr.obj$y 
     if(is.dmudr){      
           dsidr.fit<- dsidr(s=swk, q=qSum, y=y,  
                    weight=ssr.obj$weight, 
                    vmu=ssr.obj$rkpk.obj$vmu,  
                    tol=ssr.obj$call$control$tol, 
                    job=ssr.obj$call$control$job,  
                    limnla=ssr.obj$call$limnla) 
              } 
     else dsidr.fit<- ssr.obj$rkpk.obj 
## pstd=TRUE    
     if(pstd){ 
         b<- dsidr.fit$varht/(10^dsidr.fit$nlaht) 
# call dsms and dcrdr    
        if(nnull>0) dsmsV<- matrix(dsms(dsidr.fit)$sms, ncol=nnull) 
        cr<- dr <-0 
 
        for(i in 1:length(rwk)){ 
           dcrdrV<- dcrdr(dsidr.fit, rwk[[i]]) 
           cr<- cr + dcrdrV$cr %*% t(as.matrix(terms2[, i])) 
           if(nnull>0) dr<- dr + dcrdrV$dr %*% t(as.matrix(terms2[, i])) 
        } 
# collect STDs 
#  variance for the  main effect part 
      if(nnull>0){ 
        ciMain<- apply(phi, 1,  
        function(x, terms1, w) diag(terms1 %*%  
            (x*w)%*% (x*t(terms1))), terms1=terms1, w=dsmsV) 
        ciMain<- matrix(as.vector(ciMain), nrow=nrow(terms1), byrow=FALSE) 
      } 
      else ciMain<- 0 
#  variance for the rk part 
        ciRK<-  matrix( as.vector(cr) *as.vector(xi) , ncol=nobs, byrow=TRUE) 
        ciRK<- matrix(apply(ciRK, 1, sum), nrow=eval.len, byrow=FALSE) 
        Rsum<- apply(Rsum, 2, function(x, r) x[seq(1, length(x), length=r)], eval.len) 
        ciRK<- Rsum-ciRK 
# covariance between both parts 
      if(nnull>0){ 
        ciCross<-  matrix(as.vector(t(phi.new)) * as.vector(dr), ncol=nnull, byrow=TRUE) 
        ciCross<-  matrix(apply(ciCross, 1, sum), nrow=eval.len, byrow=FALSE) 
      } 
      else ciCross<- 0 
# overall covariance 
        ci<-t(ciMain) + ciRK -2* ciCross 
        ci<-ci * (ci>0) 
    } 
# caculate fits for $f$ 
    ccc<- dsidr.fit$c 
    if(nnull>0){ 
        fit<- phi%*% matDiag.prod(t(terms1), dsidr.fit$d) + apply(xi2, 2,  
               function(x, w, r) matVec.prod(matrix(x, ncol=r, byrow=TRUE), as.vector(w), left=FALSE), 
               w=ccc, r=nobs)            
    } 
    else{ 
        fit<-apply(xi2, 2,  
               function(x, w, r) matVec.prod(matrix(x, ncol=r, byrow=TRUE), as.vector(w), left=FALSE), 
               w=ccc, r=nobs)             
   } 
    
   if(ncol(fit)==1) fit<- as.vector(fit) 
    if(pstd){ 
        if(ncol(ci)==1) ci<- as.vector(ci) 
         resul[[numF]]<- list(fit=fit, pstd=sqrt(ci*b))         
    } 
    else resul[[numF]]<-list(fit.f=fit, pstd=NULL) 
  } 
  class(resul)<- c("bCI", "listbCI")
  names(resul)<- funcName 
  resul 
}     


### SNR-part
snr<-
function (formula, func, params, data = sys.parent(), start, 
    spar = "v", verbose = FALSE, control = list(), correlation = NULL, 
    weights = NULL) 
{
    thisCall <- match.call()
    y <- eval(getResponseFormula(formula)[[2]], data)
    nobs <- length(y)
    f.info <- getFunInfo(func)
    funcName <- f.info$fName
    funcArg <- f.info$f.arg
    thisCall$rk <- f.info$fRkForm
    thisCall$formF <- f.info$fNullForm
    lengthF <- length(funcName)
    if (missing(start)) 
        stop("starting values are missing")
    start.fixed <- eval(thisCall$start$params, data)
    if (verbose) {
        cat("fixed values provided: \n")
        cat(format(start.fixed))
        cat("\n")
    }
    paraModelInfo <- getParaModelMatrix(params, data, nobs)
    matrix.para <- paraModelInfo$matrix.para
    length.para <- paraModelInfo$length.para
    para.name <- paraModelInfo$para.name
    temp.index <- 1
    para.v <- NULL
    for (i in length.para) {
        temp.matrix <- matrix(matrix.para[, (temp.index):(temp.index + 
            i - 1)], ncol = i)
        para.v <- cbind(para.v, temp.matrix %*% start.fixed[(temp.index):(temp.index + 
            i - 1)])
        temp.index <- temp.index + i
    }
    para.v <- data.frame(para.v)
    names(para.v) <- para.name
    spline.f <- spline.f0 <- spline.fD <- list(lengthF)
    for (numF in 1:lengthF) {
        spline.f[[numF]] <- function(...) {
            thisFunName <- as.character(match.call()[[1]])
            thisFunOrder <- match(thisFunName, funcName)
            tau <- list(...)
            tau <- data.frame(tau)
            names(tau) <- names(x.star[[thisFunOrder]])
            result <- rkEval(thisCall$rk[[thisFunOrder]], x.star[[thisFunOrder]], 
                tau)
            result <- sumList(list.prod(result, theta.rk[[thisFunOrder]]))
            result <- matVec.prod(result, ccc[[thisFunOrder]], 
                TRUE)
            if (length(thisCall$formF) > 0 && !is.null(thisCall$formF[[thisFunOrder]])) 
                result <- result + as.vector(model.matrix2(thisCall$formF[[thisFunOrder]], 
                  tau) %*% ddd[[thisFunOrder]])
            result
        }
        spline.f0[[numF]] <- function(...) {
            thisFunName <- as.character(match.call()[[1]])
            thisFunOrder <- match(thisFunName, funcName)
            res <- data.frame(cbind(...))
            tmp.delta <- delta
            tmp.delta[[thisFunOrder]] <- res
            assign("delta", tmp.delta, envir = environment(eval(parse(text = thisFunName))))
            res
        }
        spline.fD[[numF]] <- function(...) {
            thisFunName <- as.character(match.call()[[1]])
            thisFunOrder <- match(thisFunName, funcName)
            fEst[[thisFunOrder]]
        }
    }
    if (missing(control)) 
        control <- snr.control()
    else {
        defCtrl <- snr.control()
        for (nam in names(control)) {
            if (nam == "rkpk.control") {
                for (nam.rk in names(control$rkpk.control)) {
                  defCtrl$rkpk.control[[nam.rk]] <- control$rkpk.control[[nam.rk]]
                }
            }
            else if (nam == "nls.control") {
                for (nam.nls in names(control$nls.control)) {
                  defCtrl$nls.control[[nam.nls]] <- control$nls.control[[nam.nls]]
                }
            }
            else defCtrl[[nam]] <- control[[nam]]
        }
        control <- defCtrl
    }

    convergMethod <- control$converg
    repInd <- rep(1:lengthF, rep(control$backfit, lengthF))
    f.old <- matrix(rep(1, lengthF), nrow = 1)
    coefPara.old <- c(start.fixed)
    prss.r.old <- rss <- 1
    coefF.old <- rep(mean(y), length(y))
    k <- 0
    weight <- NULL
    if (is.null(thisCall$start$f)) {
        fhat0 <- list()
        for (len in 1:lengthF) fhat0[[len]] <- 1
    }
    else fhat0 <- eval(thisCall$start$f, data)
    fhat <- fhat0 <- lapply(fhat0, function(x, y) if (length(x) > 
        1) 
        x
    else rep(x, y), y = nobs)
    rkFit <- wList <- yList <- ccc <- ddd <- theta.rk <- list()
    V <- NULL
    delta <- list()
    .env <- .GlobalEnv
    oldF <- list()
    for (i in 1:lengthF) {
        if (exists(funcName[i], envir = .env)) 
            oldF[[i]] <- eval(parse(text = funcName[i]), envir = .env)
        else oldF[[i]] <- FALSE
    }
    repeat {
        k <- k + 1
        if (verbose) 
            cat("\nIteration ", k, "\n")
        dataAug <- data.frame(data, para.v)
        for (i in 1:lengthF) {
            assign(funcName[i], spline.f0[[i]])
        }
        eval(formula[[3]], dataAug)
        dataInit <- list(lengthF)
        for (i in 1:lengthF) {
            dataInit[[i]] <- data.frame(delta[[i]])
            names(dataInit[[i]]) <- funcArg[[i]]
        }
        Q <- grid.Q <- list()
        if (k > 1) 
            grid.S <- grid.S.old
        S <- list()
        S[[lengthF + 1]] <- 1
        for (i in 1:lengthF) {
            if (length(thisCall$formF) > 0 && !is.null(thisCall$formF[[i]])) 
                S[[i]] <- model.matrix2(thisCall$formF[[i]], 
                  dataInit[[i]])
            Q[[i]] <- eval(thisCall$rk[[i]], dataInit[[i]])
            if (k > 1) 
                grid.Q[[i]] <- rkEval(thisCall$rk[[i]], dataInit[[i]], 
                  x.star[[i]])
        }
        grid.S.old <- S
        dataAug <- data.frame(data, para.v)
        for (bF in repInd) {
            fEst <- lapply(fhat, function(x) cbind(x, x, x))
            fEst[[bF]] <- cbind(fhat[[bF]], fhat[[bF]] + control$incDelta, 
                fhat[[bF]] - control$incDelta)
            for (i in 1:lengthF) assign(funcName[i], spline.fD[[i]])
            fstD <- eval(formula[[3]], dataAug)
            Dmat <- (fstD[, 2] - fstD[, 3])/control$incDelta/2
            if (bF == 1) 
                rssBF <- abs(sum((y - fstD[, 1]) * Dmat))/sqrt(sum((y - 
                  fstD[, 1])^2) * sum(Dmat^2))
            else rssBF <- max(rssBF, abs(sum((y - fstD[, 1]) * 
                Dmat))/sqrt(sum((y - fstD[, 1])^2) * sum(Dmat^2)))
            if (control$method == "GN") 
                Emat <- 0
            else {
                sndD <- (fstD[, 2] + fstD[, 3] - 2 * fstD[, 1])/control$incDelta/control$incDelta
                Emat <- sndD * (y - fstD[, 1])
            }
            if (is.null(V)) {
                uVec <- Dmat * (fstD[, 1] - y)
                DeltaMat <- abs(Dmat^2 - Emat)
                wList[[bF]] <- diag(sqrt(DeltaMat))
                yTilde <- fhat[[bF]] - uVec/DeltaMat
            }
            else {
                uVec <- Dmat * (V %*% (fstD[, 1] - y))
                DeltaMat <- Dmat * prodDiag(V, Dmat) - Emat
                wList[[bF]] <- chol.new(DeltaMat)
                yTilde <- fhat[[bF]] - solve(DeltaMat) %*% uVec
            }
            if (!is.list(Q[[bF]])) 
                Q[[bF]] <- list(Q[[bF]])
            if (length(Q[[bF]]) == 1) {
                rkFit[[bF]] <- dsidr(y = yTilde, q = Q[[bF]][[1]], 
                  s = S[[bF]], weight = wList[[bF]], tol = control$rkpk.control$tol, 
                  vmu = spar, job = control$rkpk.control$job, limnla = control$rkpk.control$limnla)
            }
            else {
                rkFit[[bF]] <- dmudr(y = yTilde, q = Q[[bF]], 
                  s = S[[bF]], weight = wList[[bF]], tol = control$rkpk.control$tol, 
                  vmu = spar, prec = control$rkpk.control$prec, maxit = control$rkpk.control$maxit, 
                  init = control$rkpk.control$init)
            }
            fhat[[bF]] <- rkFit[[bF]]$fit
            yList[[bF]] <- yTilde
            fhat0 <- fhat
        }
        for (i in 1:lengthF) {
            ccc[[i]] <- rkFit[[i]]$c
            ddd[[i]] <- rkFit[[i]]$d
            if (rkFit[[i]]$nq > 1) 
                theta.rk[[i]] <- 10^(rkFit[[i]]$theta)
            else theta.rk[[i]] <- 1
        }
        if (verbose) 
            cat("Estimate d: ", format(unlist(ddd)), "\n")
        if (k == 1) {
            grid.S <- grid.S.old
            grid.Q <- Q
        }
        f.new <- NULL
        penalty <- 0
        for (i in 1:lengthF) {
            if (is.list(grid.Q[[i]])) 
                grid.Q[[i]] <- sumList(list.prod(grid.Q[[i]], 
                  theta.rk[[i]]))
            if (is.null(grid.S[[i]])) 
                f.new <- cbind(f.new, matVec.prod(grid.Q[[i]], 
                  ccc[[i]], TRUE))
            else f.new <- cbind(f.new, matVec.prod(grid.S[[i]], 
                ddd[[i]], FALSE) + matVec.prod(grid.Q[[i]], ccc[[i]], 
                TRUE))
            if (rkFit[[i]]$nq > 1) 
                penalty <- penalty + sum(matVec.prod(grid.Q[[i]], 
                  ccc[[i]], TRUE) * ccc[[i]])
            else penalty <- penalty + sum(matVec.prod(grid.Q[[i]], 
                ccc[[i]], TRUE) * ccc[[i]]) * (10^rkFit[[i]]$nlaht)
        }
        x.star <- dataInit
        for (i in 1:lengthF) {
            assign(funcName[i], spline.f[[i]], envir = .env)
        }
        gnls.obj <- gnls(formula, data = data, params = params, 
            start = start.fixed, correlation = correlation, weights = weights, 
            control = control$nls.control)
        if (verbose) {
            cat("Estimate Coefficients: \n")
            print(gnls.obj$coef)
        }
        start.fixed <- as.vector(gnls.obj$coef)
        temp.index <- 1
        para.v <- NULL
        for (i in length.para) {
            temp.matrix <- matrix(matrix.para[, (temp.index):(temp.index + 
                i - 1)], ncol = i)
            para.v <- cbind(para.v, temp.matrix %*% start.fixed[(temp.index):(temp.index + 
                i - 1)])
            temp.index <- temp.index + i
        }
        para.v <- data.frame(para.v)
        names(para.v) <- para.name
        V <- NULL
        if (!is.null(gnls.obj$modelStruct$corStruct)) {
            V <- bdiag(corMatrix(gnls.obj$modelStruct$corStruct))
        }
        if (!is.null(gnls.obj$modelStruct$varStruct)) {
            Lambda.wei <- 1/varWeights(gnls.obj$modelStruct$varStruct)
            if (is.null(V)) 
                V <- diag(Lambda.wei^2)
            else V <- diag(Lambda.wei) %*% V %*% diag(Lambda.wei)
        }
        if (!is.null(V)) 
            V <- solve(V)
        if (convergMethod == "COEF") {
            coefPara.new <- as.vector(gnls.obj$coef)
            coefF.new <- lapply(as.list(1:lengthF), function(x, 
                m1, m2) abs(m1[, x] - m2[, x])/(1 + sum(abs(m2[, 
                x]))), m1 = f.new, m2 = f.old)
            coefF.new <- max(unlist(coefF.new))
            f.old <- f.new
            rss <- sqrt(max(coefF.new^2, (sum(abs(coefPara.new - 
                coefPara.old))/(1 + sum(abs(coefPara.old))))^2))
            coefPara.old <- coefPara.new
        }
        if (convergMethod == "PRSS") {
            prss.r.new <- sum((gnls.obj$residuals)^2) + penalty
            rss <- abs(prss.r.old - prss.r.new)/(1 + prss.r.old)
            prss.r.old <- prss.r.new
        }
        if (verbose == TRUE) 
            cat("\nConvergence Criterion:  ", rss, "\n")
        if ((rss < control$prec.out) || (k > control$maxit.out)) 
            break
    }
    for (i in 1:lengthF) {
        if (is.atomic(oldF[[i]]) && (oldF[[i]] == "FALSE")) 
            do.call("rm", list(list = funcName[i], envir = .env))
        else assign(funcName[i], oldF[[i]], envir = .env)
    }
    if (k > control$maxit.out) 
        warning("convergence not researched")
    if (!is.list(Q)) 
        Q <- list(Q)
    df.spline <- sum(unlist(lapply(rkFit, function(x) x$df)))
    dfModel <- gnls.obj$dims[["p"]] + df.spline
    sigmaMod <- gnls.obj$sigma * (nobs - gnls.obj$dims[["p"]])/(nobs - 
        dfModel)
    forCI <- list(expand.call = thisCall, data = dataInit, y = yList, 
        rkpk.obj = rkFit, s = S, q = Q, weight = wList, control = control)
    result <- list(call = match.call(), data = data, funcFitted = as.vector(f.new), 
        fitted = gnls.obj$fit, funcCoef = list(c = ccc, d = ddd), 
        coefficients = gnls.obj$coef, sigma = sigmaMod, df = list(total = nobs, 
            f = df.spline, para = gnls.obj$dims[["p"]]), forCI = forCI, 
        gnlsObj = gnls.obj, iteration = list(times = k, rss = rss))
    attr(result, "class") <- c("snr")
    result
}

print.snr<-
function(x, ...)
{
   dd <- x$gnlsObj$dims
   cat("Semi-parametric Nonlinear Regression Model Fit\n")

   cat(" Model:", deparse(as.vector(x$call$formula)), "\n")
   if(!is.null(x$call$data))
     cat(" Data:", deparse(x$call$data))
   cat("\n")
   cat(" Log-likelihood: ", format(x$gnlsObj$logLik), "\n", sep = "")   
   cat("\nCoefficients:\n")
   print(coef(x$gnlsObj))

   cat("\n")
   if(length(x$gnlsObj$modelStruct) > 0) {
		print(summary(x$gnlsObj$modelStruct))
	}    
  lambda<- NULL
  for(i in 1:length(x$forCI$rkpk.obj)){
     if(is.null(x$forCI$rkpk.obj[[i]]$theta)) lambda<- c(lambda,10^x$forCI$rkpk.obj[[i]]$nlaht)
     else lambda<- c(lambda, 10^(x$forCI$rkpk.obj[[i]]$nlaht-x$forCI$rkpk.obj[[i]]$theta))
  }
  cat("Smoothing spline:\n")  
  switch(x$forCI$rkpk.obj[[1]]$vmu,
          "v"= cat(" GCV"),
          "m"= cat(" GML"),
          "u"= cat(" UBR"))
  cat(" estimate(s) of smoothing parameter(s):", format(lambda/length(x$forCI$y[[1]])), "\n")

  cat(" Equivalent Degrees of Freedom (DF): ", format(x$df$f), "\n")
  cat("\nResidual standard error:", format(x$sigma), "\n")
  cat("Number of Observations:", length(x$forCI$y[[1]]), "\n") 
  if(!is.null(dd$groups)){
  cat("\nNumber of Groups: ")
  Ngrps <- dd$ngrps[1:dd$Q]
  if((lNgrps <- length(Ngrps)) == 1) {
# single nesting
		cat(Ngrps, "\n")
	}
	else {
# multiple nesting
		sNgrps <- 1:lNgrps
		aux <- rep(names(Ngrps), sNgrps)
		aux <- split(aux, array(rep(sNgrps, lNgrps), c(lNgrps, lNgrps))[
			!lower.tri(diag(lNgrps))])
		names(Ngrps) <- unlist(lapply(aux, paste, collapse = " %in% "))
		cat("\n")
		print(rev(Ngrps))
	}
  }
  if(x$iter[[1]] >= x$forCI$control$maxit.out) cat("\nNot Converged after", x$iter[[1]], "iterations\n")
  else cat("\nConverged after", x$iter[[1]], "iterations\n")
}

print.summary.snr<- function(x, ...){
   dd <- x$gnls.fit$dims
   cat("Semi-parametric Nonlinear Regression Model Fit by ")
   if(x$GN) cat("Gauss-Newton Method\n")
   else cat("Newton-Raphson Method\n")
   
   cat("Model:", deparse(as.vector(x$call$formula)), "\n")
   cat("Data:", deparse(x$call$data), "\n")
   cat("\n")

## summary of gnls.fit  
    print(data.frame(AIC = x$gnls.fit$AIC, BIC = x$gnls.fit$BIC, 
          logLik = x$gnls.fit$logLik, row.names = " "))   
	
    if(length(x$gnls.fit$modelStruct)) {
         cat("\n")
	 print(summary(x$gnls.fit$modelStruct))
       }
    cat("\nCoefficients:\n")
    xtTab <- as.data.frame(x$gnls.fit$tTable)
    wchPval <- match("p-value", names(xtTab))
    for(i in names(xtTab)[ - wchPval]) {
   	xtTab[, i] <- format(zapsmall(xtTab[, i]))
    }
    xtTab[, wchPval] <- format(round(xtTab[, wchPval], 4))
    if(any(wchLv <- (as.double(levels(xtTab[, wchPval])) == 0))) {
         levels(xtTab[, wchPval])[wchLv] <- "<.0001"
    }
    row.names(xtTab) <- dimnames(x$gnls.fit$tTable)[[1]]
    print(xtTab)
    if(nrow(x$gnls.fit$tTable) > 1) {
  	corr <- x$gnls.fit$corBeta
	class(corr) <- "correlation"
	print(corr, title = "\n Correlation:")
    }
    cat("\nStandardized residuals:\n")
    print(x$gnls.fit$residuals)
    cat("\n") 
  
   switch(x$vmu,
          "v"= cat("GCV"),
          "m"= cat("GML"),
          "u"= cat("UBR"))
  cat(" estimate(s) of smoothing spline parameter(s):", format(x$lamb), "\n")  
  cat("Equivalent Degrees of Freedom (DF) for spline function: ", format(x$df), "\n")
  if(is.null(x$GN))
      cat("Degrees of freedom:", dd[["N"]], "total;", format(dd[["N"]] - dd[["p"]]-x$df), "residual\n")
  cat("Residual standard error:", format(x$sigma), "\n")
  if(x$times> x$iter[[1]])
       cat("\nConverged after", x$iter[[1]], "iterations\n")
  else cat("\n Not converged after", x$iter[[1]], "iterations\n")
}
summary.snr<- function(object, ...){
## summary for snr fit   ##
  resul<-list(call=object$call,
          gnls.fit=summary(object$gnlsObj),
          times=object$forCI$control$maxit.out, 
          iter=object$iteration)
  resul$vmu<- object$forCI$rkpk.obj[[1]]$vmu
   spar <- NULL
    for (i in 1:length(object$forCI$rkpk.obj)) {
        if (object$forCI$rkpk.obj[[i]]$nq > 1) 
            spar <- c(spar, 10^(object$forCI$rkpk.obj[[i]]$nlaht-object$forCI$rkpk.obj[[i]]$theta))
        else spar <- c(spar, 10^object$forCI$rkpk.obj[[i]]$nlaht)
    }
  resul$lamb<- spar/length(object$forCI$y)
  resul$GN <- object$forCI$control$method == "GN"
  resul$df<- object$df$f
  resul$sigma<- object$sigma 
  class(resul)<- "summary.snr"
  resul
}
snr.control<-function(rkpk.control=list(job = -1, tol = 0, init = 0, limnla=c(-10, 0),
               varht=NULL, theta= NULL, prec = 1e-006, maxit = 30),            
              nls.control=list(returnObject=TRUE, maxIter=5), incDelta=0.001,
              prec.out=0.001, maxit.out=30, converg= "COEF", method="GN", backfit=5)
{
 ctrlvals.rk<-list(job = -1, tol = 0, init = 0, limnla=c(-10, 0), 
              varht=NULL, theta= NULL, prec = 1e-006, maxit = 30)
 if(!(missing(rkpk.control) || is.null(rkpk.control))){
      for(nam in names(rkpk.control)) ctrlvals.rk[[nam]]<-rkpk.control[[nam]]
 }
 ctrlvals.nls<- list(returnObject=TRUE, maxIter=5)       

 if(!(missing(nls.control) || is.null(nls.control))){
      for(nam in names(nls.control)) ctrlvals.nls[[nam]]<-nls.control[[nam]]
 }
 if(missing(prec.out) || is.null(prec.out)) prec.out<-0.0001
 if(missing(maxit.out) || is.null(maxit.out)) maxit.out<-30
 if(missing(converg) || is.null(converg)) converg<-"COEF"

 list(rkpk.control=ctrlvals.rk, nls.control=ctrlvals.nls, prec.out=prec.out,
     maxit.out=maxit.out, converg=converg, method=method, backfit=backfit, incDelta=incDelta)
}

predict.snr<- function(object, newdata=NULL, ...)
{
## Calculate fit and CI for $f$
   Call<- match.call()
   if(!inherits(object, "snr")) stop("Input doesn't match")

   ssr.obj<- object$forCI

   if(length(ssr.obj$q)>1) is.dmudr<- TRUE
   else is.dmudr<- FALSE
   if(is.null(theta<- ssr.obj$rkpk.obj$theta)) theta<- 1
   else theta<- 10^theta 

   nobs<- nrow(ssr.obj$s)
   nnull<- ssr.obj$rkpk.obj$nnull
   if(is.null(newdata)) eval.len<- nobs
   else   eval.len<- nrow(newdata)

## construct f's
   f.info<-getFunInfo(eval(object$call$func))
   funcName<- f.info$fName
   funcArg<- f.info$f.arg

   lengthF<- length(funcName)
   spline.f0<- spline.fD<- list(length(funcName))
   for(numF in 1:lengthF){            
    spline.f0[[numF]]<-function(...) {
      thisFunName<- as.character(match.call()[[1]])
      thisFunOrder<-match(thisFunName, funcName)
      res<-cbind(...)
      res.n<- nrow(res)
      res<- cbind(matrix(rep(0, (1+2*(thisFunOrder-1))*res.n), nrow=res.n), rep(1, res.n), res)
      if(length(funcName)== thisFunOrder) res 
      else
      cbind(res, matrix(rep(0, res.n*2*(length(funcName)-thisFunOrder)), nrow=res.n)) 
     }


     spline.fD[[numF]]<-function(...){
       thisFunName<- as.character(match.call()[[1]])
       thisFunOrder<-match(thisFunName, funcName)
       fEst[[thisFunOrder]]
     }
  }   
   if(!missing(newdata)){
## evaluation on newdata
     paraModel<- getParaModelMatrix(eval(ssr.obj$expand.call$params), newdata, nrow(newdata))
     start.fixed <- as.vector(object$coef)
     temp.index<-1
     para.v<-NULL     
     for(i in paraModel$length.para){
        temp.matrix<-matrix(paraModel$matrix.para[, (temp.index):(temp.index+i-1)], ncol=i)
        para.v<-cbind(para.v, 
        temp.matrix %*% start.fixed[(temp.index):(temp.index+i-1)])
        temp.index<- temp.index+i
      }
     para.v<- data.frame(para.v)  
     names(para.v)<- paraModel$para.name

     dataAug<- data.frame(newdata, para.v)

     for(i in 1:lengthF){
        assign(funcName[i], spline.f0[[i]])
     }

     delta<- eval(ssr.obj$expand.call$formula[[3]], dataAug)
     delta0<- delta[,1]

     delta<- apply(delta[,-1], 2, function(x, x0) x-x0, x0=delta0)
     dataInit<- list(lengthF)
     delta1<- NULL 
   
     fArgLen<- c(0, cumsum(unlist(lapply(funcArg, length))+1))
     for(i in 1:lengthF){
        dataTmp<- delta[, (fArgLen[i]+1):(fArgLen[i+1])]
        delta1<- cbind(delta1, dataTmp[,1])
        dataTmp<- dataTmp[,-1]/dataTmp[,1]
        dataInit[[i]]<- data.frame(dataTmp)
        names(dataInit[[i]])<- funcArg[[i]]
     }
  
     fEst<- list()
     for(i in 1:lengthF){
          if(length(ssr.obj$expand.call$formF)>0 && !is.null(ssr.obj$expand.call$formF[[i]])){
          	tmp.s<- model.matrix2(ssr.obj$expand.call$formF[[i]], dataInit[[i]]) 
#          	tmp.s<- fit.y+as.vector(tmp.s%*%ssr.obj$rkpk.obj[[i]]$d)
		tmp.s<- as.vector(tmp.s%*%ssr.obj$rkpk.obj[[i]]$d)
          }
          else tmp.s<- 0
     
        if(is.null(theta<- ssr.obj$rkpk.obj[[i]]$theta)) theta<- 1
        else theta<- 10^theta 

        tmp.q<- rkEval(ssr.obj$expand.call$rk[[i]], ssr.obj$data[[i]], dataInit[[i]])
        if(!is.list(tmp.q)) tmp.q<- list(tmp.q)

        fEst[[i]]<- tmp.s+ matVec.prod(sumList(list.prod(tmp.q,  theta)), ssr.obj$rkpk.obj[[i]]$c)
	assign(funcName[[i]], spline.fD[[i]])     
     }

     eval(ssr.obj$expand.call$formula[[3]], dataAug)
  }
   else object$gnlsObj$fit              
}

intervals.snr<- function(object, level=0.95, newdata=NULL, terms=list(), pstd=TRUE, ...)
{
## prediction of response
## fit and CI for $f$--R version
  	Call<- match.call()
	if(!inherits(object, "snr")) stop("Input doesn't match")

## spline structures
   	f.info<-getFunInfo(eval(object$call$func))
   	funcName<- f.info$fName
   	funcArg<- f.info$f.arg
   	lengthF<- length(funcName)
   	resul<- delta3<- list()

## start calculation one by one
   	for(numF in 1:lengthF){
      		ssr.obj<- object$forCI
      		ssr.obj$rkpk.obj<- ssr.obj$rkpk.obj[[numF]]
      		ssr.obj$y<- ssr.obj$y[[numF]]
      		ssr.obj$weight<- ssr.obj$weight[[numF]]
      		ssr.obj$s<- ssr.obj$s[[numF]]
      		ssr.obj$q<- ssr.obj$q[[numF]]

      		if(length(ssr.obj$q)>1) is.dmudr<- TRUE
      		else is.dmudr<-FALSE
      		if(is.null(theta<- ssr.obj$rkpk.obj$theta)) theta<- 1
      		else theta<- 10^theta 

      		nobs<- length(ssr.obj$y)
      		nnull<- ssr.obj$rkpk.obj$nnull
      		if(missing(newdata)) eval.len<- nobs
 	     	else eval.len<- nrow(newdata)

      		if(missing(terms) || is.null(terms)) terms.f<- matrix(rep(1, nnull+ length(theta)), nrow=1)
      		else {
      			if(!is.list(terms) && lengthF==1)  terms.f<-terms
        		else terms.f<-terms[[funcName[numF]]]
      		}  
      		terms.f<-as.matrix(terms.f)

      		if((ncol(terms.f) != nnull+length(theta)) && (nrow(terms.f) != nnull+length(theta)))
          		stop(" the input of terms must match ") 
   
      		if(ncol(terms.f)!= nnull+length(theta)) terms.f<-t(terms.f)
      		if(nnull>0){
     			terms1<- matrix(terms.f[, 1:(ncol(terms.f)-length(theta))], nrow=nrow(terms.f))
      			terms2<- matrix(terms.f[, (ncol(terms1)+1):(ncol(terms.f))], nrow=nrow(terms.f))
      		}
      		else{
        		terms1<- NULL
        		terms2<- terms.f
      		}

## construct f's
      		swk<- ssr.obj$s
      		qwk<- ssr.obj$q

      		if(missing(newdata)){
           		phi<- ssr.obj$s
           		rwk<- Rwk<- ssr.obj$q
          	}
      		else{
           		dataInit<- list(length(ssr.obj$data))
           		for(num in 1:length(ssr.obj$data)) dataInit[[num]]<-newdata
 
           		phi<-rwk<- Rwk<- NULL

               		if(length(ssr.obj$expand.call$formF)>0 && !is.null(ssr.obj$expand.call$formF[[numF]])){
               			tmp.s<- model.matrix2(ssr.obj$expand.call$formF[[numF]], dataInit[[numF]]) 
               			phi<- cbind(phi, tmp.s)
               		}

           		Rwk<- eval(ssr.obj$expand.call$rk[[numF]], dataInit[[numF]])
           		if(!is.list(Rwk)) Rwk<- list(Rwk)

           		rwk<- rkEval(ssr.obj$expand.call$rk[[numF]], ssr.obj$data[[numF]], dataInit[[numF]])
           		if(!is.list(rwk)) rwk<- list(rwk)
      		}

## if Null Space is null 
      		phi.new<- NULL
      		if(nnull>0){
      			for(i in 1:ncol(swk)){
           			phi.new<- cbind(phi.new, as.vector(phi[,i] %*% t(as.matrix(terms1[, i]))))
      			}
      		}
      		if(is.matrix(rwk)) rwk<- list(rwk) 
      		if(is.matrix(Rwk)) Rwk<- list(Rwk)                  

      		qSum<- Rsum<- 0
      		xi<- xi2<- 0
      		for(i in 1:length(theta)){
           		rwk[[i]]<-  rwk[[i]]*theta[i]
           		xi2<- xi2+ as.vector(rwk[[i]]) %*% t(as.matrix(terms2[, i]))
           		Rsum<- Rsum+ (as.vector(Rwk[[i]])*theta[i]) %*% t(as.matrix(terms2[, i]))
           		qSum<- qSum+ qwk[[i]]*theta[i]

           		if(!is.null(ssr.obj$weight)) rwk[[i]]<-ssr.obj$weight%*%rwk[[i]]
           		xi<- xi+ as.vector(rwk[[i]]) %*% t(as.matrix(terms2[, i]))
     		}


# call dsidr
     		y<-ssr.obj$y
     		if(is.dmudr){     
           		dsidr.fit<- dsidr(s=swk, q=qSum, y=y, 
                    		weight=ssr.obj$weight,
                    		vmu=ssr.obj$control$rkpk.control$vmu, 
                    		tol=ssr.obj$control$rkpk.control$tol,
                    		job=ssr.obj$control$rkpk.control$job, 
                    		limnla=ssr.obj$control$rkpk.control$limnla)
              	}
     		else dsidr.fit<- ssr.obj$rkpk.obj

## pstd=TRUE   
     		if(pstd){
          		b<- dsidr.fit$varht/(10^dsidr.fit$nlaht)

# call dsms and dcrdr   
          		if(nnull>0) dsmsV<- matrix(dsms(dsidr.fit)$sms, ncol=nnull)
          		cr<- dr<-0
          		for(i in 1:length(rwk)){
              			dcrdrV<- dcrdr(dsidr.fit, rwk[[i]])
              			cr<- cr + dcrdrV$cr %*% t(as.matrix(terms2[, i]))
              			dr<- dr + dcrdrV$dr %*% t(as.matrix(terms2[, i]))
          		}

# collect STDs
#  variance for the rk part
 
          		ciRK<-  matrix( as.vector(cr) *as.vector(xi) , ncol=nobs, byrow=TRUE)
          		ciRK<- matrix(apply(ciRK, 1, sum), nrow=eval.len, byrow=FALSE)
          		Rsum<- apply(Rsum, 2, function(x, r) x[seq(1, length(x), length=r)], eval.len)
          		ciRK<- Rsum-ciRK
          		if(nnull>0){
#  variance for the  main effect part
          			ciMain<- apply(phi, 1, 
          				function(x, terms1, w) diag(terms1 %*% (x*w)%*% (x*t(terms1))), terms1=terms1, w=dsmsV)
          			ciMain<- matrix(as.vector(ciMain), nrow=nrow(terms.f), byrow=FALSE)
# covariance between both parts
          			ciCross<-  matrix(as.vector(t(phi.new)) * as.vector(dr), ncol=nnull, byrow=TRUE)
          			ciCross<-  matrix(apply(ciCross, 1, sum), nrow=eval.len, byrow=FALSE)

# overall covariance
          			ci<-t(ciMain) + ciRK -2* ciCross
          		}
	  		else ci<- ciRK
          		ci<-ci * (ci>0)
     		}
# caculate fits for $f$
     		ccc<- dsidr.fit$c
     		if(nnull>0){
			fit<- phi%*% matDiag.prod(t(terms1), dsidr.fit$d) + apply(xi2, 2, 
         			function(x, w, r) matVec.prod(matrix(x, ncol=r, byrow=TRUE), as.vector(w), left=FALSE),
         			w=ccc, r=nobs)
     		}
     		else{
			fit<- apply(xi2, 2, function(x, w, r) matVec.prod(matrix(x, ncol=r, byrow=TRUE), as.vector(w), left=FALSE), w=ccc, r=nobs)
     		}
      
     		if(pstd){
        		if(ncol(ci)==1) ci<- as.vector(ci)
        		resul[[numF]]<-list(fit=fit, pstd=sqrt(ci*b))
     		}
     		else resul[[numF]]<-list(fit.f=fit, pstd=NULL)
   	}
   	class(resul)<-c("bCI", "listbCI", "list")
        names(resul)<-funcName
   	resul
}    

## SLM-part
slm<- function(formula,
               rk,
               data=sys.parent(), 
               random,
               weights=NULL,
               correlation=NULL,               
               control=list(apVar=FALSE))
##fit semi-parametric linear mixed effects models
{   
  Call <- match.call()   
  m<- match.call(expand.dots=FALSE)
  m$rk<- m$random<- m$correlation  <- m$weights<- m$scale<- m$control<- NULL
  m[[1]] <- as.name("model.frame")
  m1<- eval(m, sys.parent())
  Terms <- attr(m1, "terms")
  y <- model.extract(m1, "response")
#  y<- eval(getResponseFormula(formula)[[2]], data)

## calculate S
  S<- model.matrix(Terms, m1, NULL)
  n<-nrow(S)

## calculate Q 
	Q <- eval(Call$rk, data)
	if(is.matrix(Q)) {
		Q <- list(Q)
		Call$rk <- list(Call$rk)
	}

## form 'random'
        nq<- length(Q)
        tmp.z<- lapply(Q, chol.new)
        tmp.block<- NULL
        tmp.num<- rep(0, n)
        tmp.zz<- NULL      
        
        for(i in 1:length(tmp.z)){
          tmp.zz<-cbind(tmp.zz, tmp.z[[i]])  
          tmp.num[i+1]<- ncol(tmp.z[[i]])+tmp.num[i]         
          tmp<- list(substitute(pdIdent(~tmp.zz[, (tmp.num[i]+1):(tmp.num[i+1])]-1), list(i=i)))
          tmp.block<- c(tmp.block, tmp)
          }
        tmp.block<-lapply(tmp.block, as.formula) 
        
        tmp1<- c(ONE=tmp.block)

        if(!missing(data)) slmData<- data
        else slmData<-data.frame(ONE=rep(1, n))

        slmData[names(tmp1)]<-rep(1,n)
        tmp1<- c(tmp1, random)    
        if(any(names(tmp1)=="")){
         names(tmp1)<-paste("ONE", 1:length(tmp1), sep="")
         slmData[names(tmp1)]<-rep(1,n)
        }
        slmData$tmp.num<-tmp.num
       slmData$tmp.zz<- tmp.zz
## add groups into correlation       
     if(!is.null(correlation)){
#       cor.form<- formula.corStruct(correlation)
	cor.form<- formula(correlation)
       if(!is.null(cor.form)) cor.group<- getGroupsFormula(cor.form)
       else cor.group<- NULL
       if(!is.null(cor.group)){
         cor.group<- as.character(cor.group)[[2]]
         for(grpName in names(tmp1)[nq:1]){
            cor.group<- paste(grpName,  cor.group, sep="/")
         }
         cor.form.new<- paste(as.character(getCovariateFormula(cor.form))[[2]], "|", 
                        cor.group, sep="")
         cor.form.new<- addCorrGroup(deparse(Call$cor), cor.form.new)
         correlation<- eval(parse(text=cor.form.new))
      }
      else cor.form.new<- NULL
    }

## order by group factors ##
   if(!is.null(correlation)){
       grpfac<- NULL
       if(!is.null(names(random))) grpfac<- c(grpfac, names(random))
       if(!is.null(cor.form.new)) grpfac<- unique(c(grpfac, all.vars(getGroupsFormula(cor.form))))
       if(!is.null(grpfac)){              
         if(is.null(data)){
              dataCollect<-NULL
              for(i in grpfac) dataCollect<-cbind(dataCollect, get(i))
              dataCollect<- data.frame(dataCollect)
              names(dataCollect)<- grpfac
              grpOrder<- order.rows(dataCollect, grpfac)
            } 
         else grpOrder<- order.rows(data, grpfac)    
      }
      else grpOrder<- NULL
   }
   else grpOrder<- NULL
## call lme        		
        if(missing(data)){     
             lme.obj<- lme(formula, control=control, random=tmp1,
                   correlation=correlation, weights=weights)
        }
        else{
          lme.obj<- lme(formula, control=control, random=tmp1,
                   correlation=correlation, weights=weights, data=slmData)
        }

## get back to old cor formula
        if(!is.null(correlation) && !is.null(cor.form.new)) 
           attr(lme.obj$modelStruct$corStruct, "formula")<- cor.form		
## extract lambda
        lambda<-as.vector(coef(lme.obj$modelStruct$reStruct))
        lambda<- exp(2*lambda)
        re.est<- lme.obj$modelStruct$reStruct
        re.est[[length(re.est)]]<-NULL

## extract covariance matrix
     wei<- NULL    
     if(!is.null( lme.obj$modelStruct$corStruct)){      
       wei<- corMatrix(lme.obj$modelStruct$corStruct)
       if(!is.matrix(wei)){
          if(length(wei)==1) wei<- wei[[1]]
          else{
            wei<- bdiag(wei)
            if(!is.null(grpOrder)) wei<- wei[grpOrder,grpOrder]
           }   
        }
       cor.est <- lme.obj$modelStruct$corStruct
       }
     else cor.est<- NULL

     if(!is.null(lme.obj$modelStruct$varStruct)){
      Lambda.wei<- 1/varWeights(lme.obj$modelStruct$varStruct)
      var.est <- lme.obj$modelStruct$varStruct
      if(!is.null(wei))
        wei<- matDiag.prod(matDiag.prod(wei, Lambda.wei), Lambda.wei, left=FALSE)
      else  wei<- diag(Lambda.wei^2)
     }
    else var.est<- NULL

    re.coef<- as.vector(lme.obj$coef$random$"ONE")
    fixed.coef<-  as.vector(lme.obj$coef$fixed)
    reD<- pdMatrix(lme.obj$modelStruct$reStruct)
    if(!is.null(wei)) wei<- solve(t(chol(wei)))

## calculate zDz^t 
# get groups 
    grp<- names(tmp1)[-(1:nq)]
#    grp<- data.frame( data[grp])
    grp<- data.frame(slmData[grp])
    if((length(grp)==1) && (is.list(grp[[1]]))) 
          grpLevel<- table(grp[[1]])
    else  grpLevel<- table(grp)
    zDz<- 0

    for(ngrp in 1:ncol(grp)){
     if(ncol(grp) >1) grpi<- apply(grpLevel, ngrp, sum)
     else grpi<- grpLevel     
     ranMat<- model.matrix(reStruct(random[[ngrp]]), data=data)
     zDzi<- ranMat%*%reD[[nq+ngrp]]
     ranMat<- diagComp(t(ranMat), grpi) 
     zDzi<- diagComp(t(zDzi), grpi)
     zDzi<- t(zDzi) %*% ranMat     
     zDz<- zDz+zDzi
    }
## calculate zdz^t+W
    if(is.null(wei)) wei<- zDz+diag(n)
    else wei<- wei+zDz
    wei <- (wei + t(wei))/2
    weight <- solve(t(chol(wei)))

## call dsidr/dmudr    
   lambda.ord<-lambda<- lambda[(length(lambda)-nq+1):length(lambda)]   
   tmp.Q<- 0
   if(length(lambda) > 1) {
	for(i in 1:nq){
            lambda.ord[i]<- lambda[nq+1-i]
            tmp.Q <- tmp.Q + Q[[i]] * lambda.ord[i]
          }
	rk.obj <- dsidr(y = y, q = tmp.Q, s = S, weight= weight, vmu = "m", limnla=0)
#        lambda<- n*lambda.ord
 	lambda<- 1.0/lambda.ord/n
	}
   else {
	rk.obj <- dsidr(y = y, q = Q[[1]], s = S, 
		weight = weight, vmu = "m", limnla =  - log10(lambda[1]))
        lambda<- 1.0/lambda/n
	}
    lme.obj$call$fixed<-formula
    rk.obj<- list(call=match.call(), expand.call=Call, data=data, y=y,
                  re.est=re.est, var.est=var.est, cor.est=cor.est,
                  rkpk.obj=rk.obj, weight=weight, coef=list(re.coef=re.coef, 
                  fixed.coef=fixed.coef, d=rk.obj$d, c=rk.obj$c), vmu="m", family="gaussian",
                  residuals=lme.obj$resi[,2], fit=lme.obj$fit[,2], df=rk.obj$df,
                  lme.obj=lme.obj, s=S, q=Q, lambda=lambda, scale=scale)
    class(rk.obj)<- c("slm")
    rk.obj
}

print.slm<-function(x, ...) 
{ 
## 'print' on 'slm' class 
   dd <- x$nlme.obj$dims 
   cat("Semi-parametric linear mixed-effects model fit by ") 
   cat(ifelse(x$lme.obj$method == "REML", "REML\n", "maximum likelihood\n")) 
 
   cat("  Model:", deparse(as.vector(x$call$formula)), "\n") 
   cat("  Data:", deparse(x$call$data), "\n") 
   cat("  Log-", ifelse(x$lme.obj$method == "REML", "restricted-", ""),  
		"likelihood: ", format(x$lme.obj$logLik), "\n", sep = "")    
 
  fixF <- x$call$formula 
  if(inherits(fixF, "formula") || is.call(fixF)) { 
		cat("\nFixed:", deparse(as.vector(x$call$formula)), "\n") 
	} 
	else { 
		cat("\nFixed:", deparse(lapply(fixF, function(el) 
		as.name(deparse(as.vector(el))))), "\n") 
	} 
 
  print(fixef(x$lme.obj)) 
  cat("\n") 
  tmp<- summary(x$lme.obj$modelStruct)[[1]] 
  for(i in 1:length(x$q)) tmp[[length(tmp)]]<-NULL 
   
  if(is.null(names(x$call$random))){ 
       names(tmp)<-rep("", length(tmp)) 
   } 
 
  print(summary(tmp), sigma=x$lme.obj$sigma) 
 
  if(!is.null(x$cor.est))  
        print(x$cor.est)        
  if(!is.null(x$var.est))  
     print(x$var.est)      
 
  switch(x$rkpk.obj$vmu, 
		"v" =cat("\nGCV"), 
		"m" =cat("\nGML"), 
		"u" =cat("\nUBR"), 
		"~u"=cat("\n~UBR")) 
#	if(length(x$lambda)==1)  
#              cat(" estimate(s) of smoothing parameter(s) :", format(x$lambda), "\n") 
#        else cat(" estimate(s) of smoothing parameter(s) :", format(1.0/x$lambda), "\n") 
 	cat(" estimate(s) of smoothing parameter(s) :", format(x$lambda), "\n")

	cat("Equivalent Degrees of Freedom (DF): ", format(x$df), "\n") 
	cat("Estimate of sigma: ", format(sqrt(x$rkpk.obj$varht)), "\n")	 
                
	cat("\nNumber of Observations: ") 
	if(is.list(x$x)) 
		cat(length(x$y[[1]]), "\n") 
	else cat(length(x$y), "\n") 
} 
summary.slm<-function(object, ...){
## summary for slm fit ##
  resul<- list(call=object$call, 
                lme.fit=object$lme.obj)                
  resul$vmu<- "m"
  
#  if(length(object$lambda)>1) resul$lambda<- 1.0/object$lambda 
#  else resul$lambda<- object$lambda
  resul$lambda<- object$lambda
  resul$df<- object$rkpk.obj$df
  class(resul)<- "summary.slm" 
  resul
 }

print.summary.slm<- function(x, ...){
   cat("Semi-parametric Linear Mixed Effects Model fit\n")
   cat("  Model:", deparse(as.vector(x$call$formula)), "\n")
   cat("  Data:", deparse(x$call$data), "\n")
   cat("\n")
## summary of lme.fit 
   lme.sum<- summary(x$lme.fit)
   lme.sum$call$fixed<-x$call$formula
   lme.sum$call$data<-x$call$data
   length.re<- length(lme.sum$modelStruct["reStruct"][[1]])
   for(term in 1:length(x$lambda)){
      lme.sum$modelStruct["reStruct"][[1]][[length.re-term+1]]<-NULL
      lme.sum$coef$random[[term]]<-NULL
     }

   lme.sum$dims$ngrps<- lme.sum$dims$ngrps[-c(length.re:(length.re-length(x$lambda)+1))]
   lme.sum$dims$Q<- lme.sum$dims$Q-length(x$lambda)
   print(lme.sum)
## summary of f
   cat("\nSmoothing spline:\n") 
    switch(x$vmu,
          "v"= cat(" GCV"),
          "m"= cat(" GML"),
          "u"= cat(" UBR"))

  cat(" estimate(s) of smoothing parameter(s):", format(x$lamb), "\n")  
  cat(" Equivalent Degrees of Freedom (DF): ", format(x$df), "\n")
}

predict.slm<-
function (object, newdata = NULL, ...) 
{
    if (!missing(newdata)) {
        lmeCoef <- object$lme.obj$coef
        coef.ran <- lmeCoef$random
        struct.ran <- reStruct(eval(object$call$random))
        form.ran <- formula(struct.ran)
        name.ran <- names(form.ran)
        grp.ran <- getGroupsFormula(struct.ran)
        if (is.null(object$data)) {
            grpName <- matrix(all.vars(grp.ran[[2]]), nrow = 1)
            grp.old <- data.frame(apply(grpName, 2, get))
            names(grp.old) <- grpName
        }
        else grp.old <- data.frame(getGroups(eval(object$data), 
            grp.ran))
        grp.ran <- getGroups(newdata, grp.ran)
        dMat.ran <- model.matrix(struct.ran, newdata)

## =====>this part is for R only
	mat.nam<- attr(dMat.ran, "nams")
	mat.len<- lapply(mat.nam, length)
	mat.ind<- unlist(sapply(1:length(mat.nam), function(x, y) rep(x, y[[x]]), y=mat.len))
	mat.ind<- split(1:length(mat.ind), mat.ind)
	dMat.ran<- rev(lapply(mat.ind, function(x, y) y[,x], y=dMat.ran))
## <==============
        if (inherits(grp.ran, "factor")) {
            ord.ind <- order(grp.ran)
            grp.ran <- sort(grp.ran)
            nobs.in <- list(table(grp.ran))
            nobs.gp <- length(nobs.in[[1]])
            grp.ran <- data.frame(grp.ran)
        }
        else {
            ord.ind <- do.call("order", grp.ran)
            nobs.in <- apply(grp.ran, 2, table)
            nobs.gp <- unlist(lapply(nobs.in, length))
        }
        ord.org <- as.character(1:nrow(newdata))
        ord.org <- match(ord.org, ord.org[ord.ind])

        if (is.matrix(dMat.ran)) 
            dMat.ran <- list(dMat.ran)
        coef.ran <- rev(coef.ran)
        length(coef.ran) <- length(dMat.ran)
        coef.ran <- rev(coef.ran)
        for (i in 1:length(coef.ran)) {
            existed <- match(unique(grp.ran[, i]), sort(unique(grp.old[, 
                i])))
            coef.ran[[i]] <- coef.ran[[i]][existed, ]
        }
        dMat.ran2 <- as.list(1:length(dMat.ran))
        dMat.ran2 <- lapply(dMat.ran2, function(x, mat, ind) diagComp(t(mat[[x]]), 
            ind[[x]]), mat = dMat.ran, ind = nobs.in)
        fit.ran <- as.list(1:length(dMat.ran))
        fit.ran <- lapply(fit.ran, function(x, mat, coe) t(mat[[x]]) %*% 
            as.vector(t(coe[[x]])), mat = dMat.ran2, coe = coef.ran)
        fit.ran <- matrix(unlist(fit.ran), ncol = length(fit.ran))
        fit.f <- intervals(object, newdata = newdata, pstd = FALSE)$fit
        fitted <- cbind(fit.f, fit.ran)
        fitted <- as.data.frame(t(apply(fitted, 1, cumsum)))
        dimnames(fitted)[[2]] <- c("fixed", name.ran)
    }
    else {
        fitted <- object$lme.obj$fit[, -c(1:length(object$q))]
        dimnames(fitted)[[2]][[1]] <- "fixed"
    }
    fitted
}

   
intervals.slm<-
function(object,level=0.95,  newdata = NULL, terms, pstd = TRUE, ...)
{
## calculate posterior STD and fits 
## of SS ANOVA components from slm object
	Call <- match.call()
	nobs <- length(object$y)
	if(length(theta <- object$lambda) == 1)
		theta <- 1
#	else theta <- theta/nobs
        else theta<- 1/(theta*nobs)

	nnull <- object$rkpk.obj$nnull
	if(is.null(newdata))
		eval.len <- nobs
	else eval.len <- length(newdata[[1]])
        data.used <- object$data
	if(missing(terms) || is.null(terms))
		terms <- rep(1, nnull + length(theta))
	terms <- as.matrix(terms)
	if((ncol(terms) != nnull + length(theta)) && (nrow(terms) != nnull + 
		length(theta)))
		stop(" the input of terms must match ")
	if(ncol(terms) != nnull + length(theta))
		terms <- t(terms)
	terms1 <- matrix(terms[, 1:(ncol(terms) - length(theta))], nrow = nrow(
		terms))
	terms2 <- matrix(terms[, (ncol(terms1) + 1):(ncol(terms))], nrow = nrow(
		terms))	

# calculate S and Q and r's
	swk <- object$s
	qwk <- object$q
	if(is.null(newdata)) {
		phi <- object$s
		rwk <- Rwk <- object$q
	}
	else {
		rwk <- rkEval(object$expand.call$rk, data.used, newdata)
             	
#	## if subset
#		if(!is.null(eval.data(object$call$subset, object$call$data))
#			) {
#			for(compon in 1:length(rwk)) {
#				rwk[[compon]] <- rwk[[compon]][eval.data(
#				  object$call$subset, object$call$data),  ]
#			}
#		}
#
## end-of-subset
		phi <- model.matrix(eval(object$expand.call$formula[-2]), data=newdata)
		if(nrow(phi) == 1)
			phi <- matrix(rep(as.vector(phi), dim(newdata)[1]), 
				nrow = dim(newdata)[1])
		Rwk <- eval(object$call$rk, newdata)
		if(is.matrix(Rwk))
			Rwk <- list(Rwk)
	}
	phi.new <- NULL
	for(i in 1:ncol(swk)) {
		phi.new <- cbind(phi.new, as.vector(phi[, i] %*% t(as.matrix(
			terms1[, i]))))
	}
	if(is.matrix(rwk))
		rwk <- list(rwk)
	if(is.matrix(Rwk))
		Rwk <- list(Rwk)
	qSum <- Rsum <- 0
	xi <- xi2 <- 0
	for(i in 1:length(theta)) {
		rwk[[i]] <- rwk[[i]] * theta[i]
		xi2 <- xi2 + as.vector(rwk[[i]]) %*% t(as.matrix(terms2[, i]))
		Rsum <- Rsum + (as.vector(Rwk[[i]]) * theta[i]) %*% t(as.matrix(
			terms2[, i]))
		qSum <- qSum + qwk[[i]] * theta[i]
		if(!is.null(object$weight))
			rwk[[i]] <- object$weight %*% rwk[[i]]
		xi <- xi + as.vector(rwk[[i]]) %*% t(as.matrix(terms2[, i]))
	}
# call dsidr
	y <- object$y
	if(!is.null(newdata)) {
		dsidr.fit <- dsidr(s = swk, q = qSum, y = y, weight = object$
			weight, vmu = "m")
	}
	else dsidr.fit <- object$rkpk.obj
	if(pstd) {
		b <- dsidr.fit$varht/(10^dsidr.fit$nlaht)	
	# call dsms and dcrdr   
		dsmsV <- matrix(dsms(dsidr.fit)$sms, ncol = nnull)
		cr <- 0
		dr <- 0
		for(i in 1:length(rwk)) {
			dcrdrV <- dcrdr(dsidr.fit, rwk[[i]])
			cr <- cr + dcrdrV$cr %*% t(as.matrix(terms2[, i]))
			dr <- dr + dcrdrV$dr %*% t(as.matrix(terms2[, i]))
		}
# collect STDs
#  variance for the  main effect part
		ciMain <- apply(phi, 1, function(x, terms1, w)
		diag(terms1 %*% (x * w) %*% (x * t(terms1))), terms1 = terms1, 
			w = dsmsV)
		ciMain <- matrix(as.vector(ciMain), nrow = nrow(terms), byrow
			 = FALSE)	#  variance for the rk part
		ciRK <- matrix(as.vector(cr) * as.vector(xi), ncol = nobs, 
			byrow = TRUE)
		ciRK <- matrix(apply(ciRK, 1, sum), nrow = eval.len, byrow = FALSE)
		Rsum <- apply(Rsum, 2, function(x, r)
		x[seq(1, length(x), length = r)], eval.len)
		ciRK <- Rsum - ciRK	# covariance between both parts
		ciCross <- matrix(as.vector(t(phi.new)) * as.vector(dr), ncol
			 = nnull, byrow = TRUE)
		ciCross <- matrix(apply(ciCross, 1, sum), nrow = eval.len, 
			byrow = FALSE)	# overall covariance
		ci <- t(ciMain) + ciRK - 2 * ciCross
		ci <- ci * (ci > 0)
	}
# caculate fits
	ccc <- dsidr.fit$c
	fit <- phi %*% matDiag.prod(t(terms1), dsidr.fit$d) + apply(xi2, 2, 
		function(x, w, r)
	matVec.prod(matrix(x, ncol = r, byrow = TRUE), as.vector(w), left = FALSE), w
		 = ccc, r = nobs)
	if(ncol(fit) == 1)
		fit <- as.vector(fit)
	if(pstd) {
		if(ncol(ci) == 1)
			ci <- as.vector(ci)
		resul<-list(fit = fit, pstd = sqrt(ci * b))
	}
	else resul<-list(fit = fit, pstd = NULL)
  class(resul)<-c("bCI", "list")
  resul 
}

### SNM-part
### SNM-part
snm<- function(formula,
             func,
             data=sys.parent(),
             fixed, 
             random=fixed,
             groups,
             start,
	     spar="v",
             verbose=FALSE,
             method="REML",                                      
             control=NULL,
             correlation=NULL,
             weights=NULL) 
{ 
## snm for R  
## allow multiple (linear) f's
  Call <- match.call()

## figure out response and grouping var's
  y<- eval(getResponseFormula(formula)[[2]], data)
  nobs <- length(y)
  if(missing(groups)){
       groups<-getGroupsFormula(reStruct(random))
   }
  if(missing(data)){
    grpName<-all.vars(groups[[2]])
    grp<- as.list(grpName)
    grp<-data.frame(lapply(grp, get))
    names(grp)<- grpName
    if(ncol(grp)==1) grp<-as.factor(grp[,1])
  }
  else  grp<-getGroups(data, groups)

## order data based on groups   
  if(inherits(grp, "factor")){
     ord.ind<-order(grp)
      grp<- sort(grp)
      nobs.in <- list(table(grp))
      nobs.gp <- length(nobs.in[[1]])
     }
  else{
      ord.ind<- do.call("order", grp) 
      nobs.in<- lapply(grp, function(x) table(sort(x)) )
      nobs.gp<- unlist(lapply(nobs.in, length))
     }    
   
  ord.org<- as.character(1:nobs)
  ord.org<- match(ord.org, ord.org[ord.ind])

##extract the start values for fixed effects
  if(!missing(start)) start.fixed<-eval(Call$start)
  else start.fixed<-mean(y)
#  if(missing(start)) stop("Start values are required")
#  else{
#      if(is.list(start)){
#                start.fixed<-start$fixed
#                start.random<- start$random
#      }
#      else if(is.vector(start)) start.fixed<-eval(Call$start)
#  }
  
  if(verbose==TRUE) cat("Starting values provided: ", start.fixed, "\n")

## extract effect names
## fixed effects

  fixModelInfo<- getParaModelMatrix(fixed, data, nobs)
  matrix.fix<- fixModelInfo$matrix.para
  length.fix<- fixModelInfo$length.para
  nameFixed<- fixModelInfo$para.name
# order the matrix
  matrix.fix<- matrix.fix[ord.ind, , drop=FALSE]

##random effects  
  random<- reStruct(random)
 # re.form<- (formula.reStruct(random))[[1]]
  re.form<- (formula(random))[[1]]

  nameRandom<- unlist(lapply(re.form, function(x) as.character(getResponseFormula(x)[[2]])))
  if(missing(data)){
   dataTmp<-data.frame(matrix(rep(1, nobs*length(nameRandom)), ncol=length(nameRandom)))
   names(dataTmp)<-nameRandom
   }
  else{
    dataTmp<- data
    dataTmp[nameRandom]<-rep(1, nrow(dataTmp))
  }
  matrix.random<- model.matrix(random, dataTmp)
  length.random<- unlist(lapply(re.form, 
          function(x, z) ncol(model.matrix2(getCovariateFormula(x), z)), z=dataTmp))

# order matrix.random
  matrix.random<- matrix.random[ord.ind, , drop=FALSE]

## combine fixed and random effects
  paraName<- c(nameFixed, nameRandom[match(nameRandom, nameFixed, nomatch=0)==0])

## initial values to evaluate deltas
  temp.index<-1
  para.v<-NULL

  for(i in length.fix){
     temp.matrix<-matrix(matrix.fix[, (temp.index):(temp.index+i-1)], ncol=i)
     para.v<-cbind(para.v, 
          temp.matrix %*% start.fixed[(temp.index):(temp.index+i-1)])
     temp.index<- temp.index+i
    }
  if((length(paraName)-length(nameFixed))>0)
      para.v <- cbind(para.v, 
         matrix( rep(0, (length(paraName)-length(nameFixed))*nobs), nrow=nrow(para.v)))
  para.v <- data.frame(para.v)
  names(para.v) <- paraName

## the following allows to input random effects
#  para.v<- getParaValue(length.fix, length.random, matrix.fix, 
#        matrix.random,start.fixed, start.random, paraName, nameRandom, 
#         nobs, nameFixed, nobs.in)
#print(para.v)
  start.cov <- diag(rep(1, ncol(matrix.random)))

## get info for "f"
  f.info<-getFunInfo(func)
  funcName<- f.info$fName
  funcArg<- f.info$f.arg
  Call$rk<- f.info$fRkForm
  Call$formF<- f.info$fNullForm
  lengthF<- length(funcName)

## construct spline $f$ 
  spline.f<- spline.f0<- list(lengthF)
  for(numF in 1:lengthF){
      spline.f[[numF]]<-function(...){
         thisFunName<- as.character(match.call()[[1]])
         thisFunOrder<-match(thisFunName, funcName)
         tau<- data.frame(...)
         names(tau)<- names(x.star[[thisFunOrder]])
         result<- rkEval(Call$rk[[thisFunOrder]], x.star[[thisFunOrder]], tau)	
         result<- sumList(list.prod(result, 
             theta.rk[(q.index[thisFunOrder]+1):(q.index[thisFunOrder+1])]))
         result<- matVec.prod(result, ccc*delta1[,thisFunOrder], TRUE)
         if(length(Call$formF)>0 && !is.null(Call$formF[[thisFunOrder]])) 
           result<- result+ as.vector(model.matrix2(Call$formF[[thisFunOrder]], tau)%*%
                   ddd[(d.index[thisFunOrder]+1):(d.index[thisFunOrder+1])])
         result
       }
            
     spline.f0[[numF]]<-function(...) {
         thisFunName<- as.character(match.call()[[1]])
         thisFunOrder<-match(thisFunName, funcName)
         res<-data.frame(...)        
         res<-sapply(res, as.numeric)
         res.n<- nrow(res)
         res<- cbind(matrix(rep(0, (1+2*(thisFunOrder-1))*res.n), nrow=res.n), rep(1, res.n), res)
         if(length(funcName)== thisFunOrder) res 
         else
          cbind(res, matrix(rep(0, res.n*2*(length(funcName)-thisFunOrder)), nrow=res.n)) 
        }
  }   


  fArgLen<- c(0, cumsum(unlist(lapply(funcArg, length))+1))

## match control options
  if(missing(control)) control<- snm.control()
  else control<- snm.control(control$rkpk.control, nlme.control=control$nlme.control, 
       prec.out=control$prec.out, maxit.out=control$maxit.out, converg=control$converg)

  convergMethod<- control$converg  

##   iteration starts
  V<- NULL
  coefPara.old <-c(start.fixed, as.vector(diag(start.cov)))
  prss.r.old<- rss<- 1
  coefF.old<-rep(mean(y), nobs)
  k<- 0

## handle incosistency of environment in nlme
## save variables already existing in the workpace 
## which will be restored later. This is kind of ugly,
## but not easy to go around.
#  .env<- environment(nlme)
  .env<- .GlobalEnv

  oldF<-list()
  for(i in 1:lengthF){
     if(exists(funcName[i], envir=.env)) oldF[[i]]<-eval(parse(text=funcName[i]), envir=.env)
     else oldF[[i]]<-FALSE 
  } 

  repeat{
    k<- k+1
    if(verbose==TRUE) cat("\nIteration ", k, "\n")     
    dataAug<- data.frame(data, para.v[ord.org, , drop=FALSE])

## get delta1 and delta2 ##
    for(i in 1:lengthF){
      assign(funcName[i], spline.f0[[i]], envir=.env)      
    } 
    delta<- eval(formula[[3]], dataAug)

    delta0<- delta[,1]
    delta<- apply(delta[,-1], 2, function(x, x0) x-x0, x0=delta0)
    dataInit<- list(lengthF)
    delta1<- NULL 
    for(i in 1:lengthF){
       dataTmp<- delta[, (fArgLen[i]+1):(fArgLen[i+1])]
       delta1<- cbind(delta1, dataTmp[,1])
       dataTmp<- dataTmp[,-1]/dataTmp[,1]
       dataInit[[i]]<- data.frame(dataTmp)
       names(dataInit[[i]])<- funcArg[[i]]
    }

##calculate Q and S
    Q<- grid.Q<- NULL
    if(k==1){
      d.index<-  q.index<- 0
      S<- NULL 
      for(i in 1:lengthF){   
         if(length(Call$formF)>0 && !is.null(Call$formF[[i]])) 
            S<- cbind(S, delta1[,i]*model.matrix2(Call$formF[[i]], dataInit[[i]]))
         d.index<- c(d.index, ncol(S))
         tmp.q<- eval(Call$rk[[i]], dataInit[[i]])
         if(!is.list(tmp.q)) tmp.q<- list(tmp.q)
         tmp.q<- lapply(tmp.q, function(x, x0)  
                  matDiag.prod(matDiag.prod(x, x0, TRUE), x0, FALSE), x0=delta1[,i])
         Q<- c(Q, tmp.q)
         q.index<-c(q.index, length(Q))
         }
     }
    else{
#      grid.S<- S
      grid.S<-S<- NULL
      for(i in 1:lengthF){
         if(length(Call$formF)>0 && !is.null(Call$formF[[i]])){ 
            S<- cbind(S, delta1[, i]*model.matrix2(Call$formF[[i]], dataInit[[i]]))        
            grid.S<- cbind(grid.S, model.matrix2(Call$formF[[i]], x.star[[i]]))
	    }
         tmp.q<- eval(Call$rk[[i]], dataInit[[i]])
         if(!is.list(tmp.q)) tmp.q<- list(tmp.q)        
         Q<- c(Q, lapply(tmp.q, function(x, x0)
                 matDiag.prod(matDiag.prod(x, x0, TRUE), x0, FALSE), x0=delta1[,i]))        
         tmp.q<- rkEval(Call$rk[[i]], dataInit[[i]], x.star[[i]])
         if(!is.list(tmp.q)) tmp.q<- list(tmp.q)
         grid.Q<- c(grid.Q, 
            lapply(tmp.q, function(x, x0) (diag(x0)%*%x), x0=delta1[,i]))
        }
     }
#    if(is.null(S))
#      S<- matrix(rep(1, nobs), ncol=1)
    
## rkpk step begins 
    if((!is.list(Q)) || (length(Q)==1)){
       if(length(Q)==1) Q<- Q[[1]]
       rk.obj<- dsidr(y=y-delta0, q=Q, s=S,
                   tol=control$rkpk.control$tol,  vmu=spar,
                   job=control$rkpk.control$job,  limnla=control$rkpk.control$limnla)       
       }
    else{
          rk.obj <- dmudr(y=y-delta0, q=Q, s=S,
                   tol=control$rkpk.control$tol,  vmu=spar,
                   prec=control$rkpk.control$prec,  maxit=control$rkpk.control$maxit,
                   init=control$rkpk.control$init)  
       }
    if(is.null(rk.obj$theta)) theta.rk<- 1    
    else theta.rk<- 10^(rk.obj$theta)   

    if(verbose){
       if(length(theta.rk)==1)   
           cat("smoothing parameters: ", format(10^(rk.obj$nlaht)), "\n")
       else cat("smoothing parameters: ", format(theta.rk), "\n")
     }

    ccc<- rk.obj$c
    ddd<-rk.obj$d

    if(length(d.index)==1) ddd<- 0
    if(verbose) cat("d: ", format(ddd), "\n")
#if(verbose) cat("d: ", format(rk.obj$d), "\n")
## calculate fit on grid points
    if(k>1){
       if(is.list(grid.Q)) grid.Q<- sumList(list.prod(grid.Q, theta.rk))
       fit.new<- matVec.prod(grid.Q, ccc, TRUE)
       if(!is.null(grid.S)) fit.new<- matVec.prod(grid.S, ddd, FALSE)+fit.new
    }
    else fit.new<- rk.obj$fit
        
## rkpk step ends

    if(convergMethod=="PRSS"){
      lambda<- 10^(rk.obj$nlaht)
      if(is.list(Q)) penalty<- sum(matVec.prod(sumList(list.prod(Q, theta.rk)), rk.obj$c, TRUE)*rk.obj$c)
      else penalty<- sum(matVec.prod(Q, rk.obj$c, TRUE)*rk.obj$c)*lambda
     }
 
## begin nlme step
    x.star<- dataInit 
    for(i in 1:lengthF){
         assign(funcName[i], spline.f[[i]], envir=.env)
       }
    if(k==1){
        nlme.obj<-nlme(formula, fixed=fixed,
           random=random, groups=groups,
           data=data, method=method,
           start=list(fixed=start.fixed, D=start.cov),
           correlation=correlation, weights=weights,
           control=control$nlme.control)
        }

    else{
       nlme.obj<-nlme(formula, fixed=fixed,
           random=random, groups=groups,
           data=data, method=method,
           start=list(fixed=start.fixed, random=start.random,
           D=start.cov), correlation=correlation, weights=weights,
           control=control$nlme.control)
      }

   if(verbose==TRUE){
     cat("Estimate of Coef: \n ")
     print(nlme.obj$coef)
     }

## update the starting values for next-step nlme
    coef.nlme<- nlme.obj$coef
    start.cov <- pdMatrix(nlme.obj$modelStruct$reStruct)
    if(verbose==TRUE) {
       cat("Estimate of D: \n ")
       print(start.cov)
      }
   start.fixed <- as.vector(coef.nlme$fixed)
   start.random<- nlme.obj$coef$random
## start of testing
## force random sum to zero for the first two step
## to be refined
##   if(k<=5) 
#start.random<- lapply(start.random, function(x) apply(x,2, function(y) y-c(rep(mean(y[1:11]),11), rep(mean(y[-c(1:11)]),9))  ))
#start.random<- lapply(start.random, function(x) apply(x,2, function(y) y-mean(y)  ))
## end of testing

## prepare for the evaluations of 3 delta's  
   para.v<- getParaValue(length.fix, length.random, matrix.fix, 
        matrix.random,start.fixed, start.random, paraName, nameRandom, 
         nobs, nameFixed, nobs.in)

##  Two ways to control the convergence
   if(convergMethod=="COEF"){    
       coefPara.new<-  c(start.fixed, as.vector(unlist(lapply(start.cov, diag))))
       coefF.new<-  sum(abs(fit.new-coefF.old))/(1+sum(abs(coefF.old)))
       coefF.old<- fit.new           
       rss<- max( coefF.new^2, 
             sum(abs(coefPara.new-coefPara.old))/(1+sum(abs(coefPara.old))))
      coefPara.old<- coefPara.new
     }

    else if(convergMethod=="PRSS") {
         prss.r.new<- sum(( nlme.obj$residuals[,2])^2) + penalty
         rss<- abs(prss.r.old-prss.r.new)/(1+prss.r.old)              
         prss.r.old <- prss.r.new      
      }

    if(verbose==TRUE)
      cat("\nConvergence Criterion:  ", rss, "\n")    

## If weight!=NULL ##  
   V<- NULL
   if(!is.null(nlme.obj$modelStruct$corStruct)){
      V<- bdiag(corMatrix(nlme.obj$modelStruct$corStruct))
     }

   if(!is.null(nlme.obj$modelStruct$varStruct)){
      if(is.null(V))
          V<- diag(varWeights(nlme.obj$modelStruct$varStruct))
      else{
         tmp.V<- 1/varWeights(nlme.obj$modelStruct$varStruct)
#         V<- matDiag(matDiag.prod(V, tmp.V, TRUE), tmp.V, FALSE)
         V<- matDiag.prod(matDiag.prod(V, tmp.V, TRUE), tmp.V, FALSE)
         V<- solve(t(chol( (V+t(V))/2)))
        }
     }
   else{
     if(!is.null(V)) V<- solve(t(chol( (V+t(V))/2)))      
     }
   if( (rss<control$prec.out) || (k > control$maxit.out)) break 
  }
## iteration ends

  if(k >= control$maxit.out) warnings("Convergence not reached!")
## numerical derivatives
  numLevel<- length(nobs.gp)
  parDeriv<-list(numLevel)
  incDelta<- 0.001

  paraValInc<- getParaValue(length.fix, length.random, matrix.fix, 
          matrix.random,start.fixed, start.random, paraName, nameRandom, 
          nobs, nameFixed, nobs.in)

  dataAugInc<- data.frame(data, paraValInc[ord.org, ,drop=FALSE])
  derivFit<- eval(formula[[3]], dataAugInc)
  for(numList in 1:numLevel){
     tmpParDeriv<- NULL
     startRandInc<- start.random
     for(i in 1:ncol(start.random[[numList]])){       
       startRandInc[[numList]][,i]<-start.random[[numList]][,i]+incDelta

       paraValInc<- getParaValue(length.fix, length.random, matrix.fix, 
          matrix.random,start.fixed, startRandInc, paraName, nameRandom, 
          nobs, nameFixed, nobs.in)
       dataAugInc<- data.frame(data, paraValInc[ord.org, , drop=FALSE])
       deltaInc<- eval(formula[[3]], dataAugInc)

       tmpParDeriv<- cbind(tmpParDeriv, (deltaInc-derivFit)/incDelta)
      startRandInc[[numList]]<-start.random[[numList]]
     }
     parDeriv[[numList]]<-tmpParDeriv
  }

## calculate y.star=y+ Z*b0
  Z.new<- as.list(1:numLevel)
  Z.new<- lapply(Z.new, 
         function(x, mat, nGrp) diagComp(t(mat[[x]]), nGrp[[x]]), mat=parDeriv, nGrp=nobs.in)
  y.star<- y+apply(matrix(1:length(start.random), ncol=1), 2, 
     function(x, mat, coefRand) t(mat[[x]])%*%as.vector(t(coefRand[[x]])), mat=Z.new, coefRand=start.random)

  D<- V<- as.list(1:numLevel)
  D<- lapply(D, function(x, mat,ind) 
  diagComp(matrix(rep(mat[[x]], ind[x]), nrow=nrow(mat[[x]])), rep(nrow(mat[[x]]), ind[[x]])),
          mat=start.cov, ind=nobs.gp)  
  V<-lapply(V, function(x, mat1, mat2) t(mat1[[x]])%*%mat2[[x]]%*%mat1[[x]], mat1=Z.new, mat2=D)
  V<- sumList(V)
## correltation structure

  if(!is.null(nlme.obj$modelStruct$corStruct))
   Lamda.cor<- bdiag(corMatrix(nlme.obj$modelStruct$corStruct))
#   Lamda.cor<- diag.matrix(corMatrix(nlme.obj$modelStruct$corStruct))
  else Lamda.cor<- diag(rep(1, ncol(V)))

## varFunction structure
  if(!is.null(nlme.obj$modelStruct$varStruct)){
      Lambda.wei<- 1/varWeights(nlme.obj$modelStruct$varStruct)
      Lamda.cor<- diag(Lambda.wei) %*% Lamda.cor %*% diag(Lambda.wei)
     }

  V<- V+ Lamda.cor
  V<- (V+t(V))/2
  V<- solve(t(chol(V)))

## rkpk step
  if(!is.list(Q)) Q<- list(Q)
  if(length(Q)==1){
       rk.obj2<- dsidr(y=y.star, q=Q[[1]], s=S, weight=V,
             limnla=control$rkpk.control$limnla, vmu=spar, 
             job=control$rkpk.control$job, tol=control$rkpk.control$tol)  
     }
  else{
       rk.obj2 <- dmudr(y=y.star, q=Q, s=S, weight=V, 
                   tol=control$rkpk.control$tol, maxit=control$rkpk.control$maxit, 
                   vmu=spar, prec=control$rkpk.control$prec, init=control$rkpk.control$init)
      } 
   
  if(is.null(rk.obj2$nq)) rk.obj2$nq<-1
  forCI<-list(expand.call=Call, delta1=delta1, data=dataInit, control=control, 
                  weight=V, y=y.star, rkpk.obj=rk.obj2, q=Q, s=S, grpfac=grp)

## restore variables covered
 for(i in 1:lengthF){
     if(is.atomic(oldF[[i]]) && (oldF[[i]]=="FALSE")) 
        do.call("rm", list(list=funcName[i], envir=.env)) 
     else  assign(funcName[i], oldF[[i]], envir=.env)                    
 }  

## adjusted sigma
 if(method=="REML") sigmaDF<-length(nlme.obj$coef$fixed)
 else sigmaDF<-0
 sigma<- nlme.obj$sigma*sqrt(nobs-sigmaDF)/sqrt(nobs-sigmaDF-rk.obj$df)
 nlme.obj$sigma<-sigma
  result<-list(call=match.call(),
              data=data, 
              vmu=spar,
              lambda=rk.obj$nlaht,
              theta=rk.obj$theta,
              funcFitted=as.vector(fit.new),
              fitted=nlme.obj$fit[,2],
              funcCoef=list(c=ccc, d=ddd), 
              coefficients=nlme.obj$coef,
              funcDf=rk.obj$df,
              sigma=sigma,
              nlmeObj=nlme.obj,  
              iteration=list(times=k, converg=rss), 
              forCI=forCI)
  attr(result, "class")<- c("snm")
  result
}


intervals.snm<-function (object, level = 0.95, newdata = NULL, terms, pstd = TRUE, 
    ...) 
{
    Call <- match.call()
    if (!inherits(object, "snm")) 
        stop("Input doesn't match")
    ssr.obj <- object$forCI
    if (length(ssr.obj$q) > 1) 
        is.dmudr <- TRUE
    else is.dmudr <- FALSE
    if (is.null(theta <- object$theta)) 
        theta <- 1
    else theta <- 10^theta
    nobs <- ssr.obj$rkpk.obj$nobs
    nnull <- ssr.obj$rkpk.obj$nnull
    if (is.null(newdata)) 
        eval.len <- nobs
    else eval.len <- nrow(newdata)
    if (missing(terms) || is.null(terms)) 
        terms <- rep(1, nnull + length(theta))
    terms <- as.matrix(terms)
    if ((ncol(terms) != nnull + length(theta)) && (nrow(terms) != 
        nnull + length(theta))) 
        stop("the input of terms must match")
    if (ncol(terms) != nnull + length(theta)) 
        terms <- t(terms)
    if (nnull > 0) {
        terms1 <- matrix(terms[, 1:(ncol(terms) - length(theta))], 
            nrow = nrow(terms))
        terms2 <- matrix(terms[, (ncol(terms1) + 1):(ncol(terms))], 
            nrow = nrow(terms))
    }
    else {
        terms1 <- NULL
        terms2 <- terms
    }
    f.info <- getFunInfo(eval(object$call$func))
    funcName <- f.info$fName
    funcArg <- f.info$f.arg
    lengthF <- length(funcName)
    swk <- ssr.obj$s
    qwk <- ssr.obj$q

    if(missing(newdata))  dataInit<- ssr.obj$data
    else {
        dataInit <- list(length(ssr.obj$data))
        for (num in 1:length(dataInit)) dataInit[[num]] <- newdata
        }
    
     phi <- rwk <- Rwk <- NULL
     oldCoefD <- ssr.obj$rkpk.obj$d
     thetaUsed <- theta
     for (i in 1:lengthF) {
            if (length(ssr.obj$expand.call$formF) > 0 && !is.null(ssr.obj$expand.call$formF[[i]])) {
                tmp.s <- model.matrix2(ssr.obj$expand.call$formF[[i]], 
                  data = dataInit[[i]])
                phi <- cbind(phi, tmp.s)
                oldCoefD <- oldCoefD[-c(1:ncol(tmp.s))]
            }
     	    tmp.q <- eval(ssr.obj$expand.call$rk[[i]], dataInit[[i]])
            if (!is.list(tmp.q)) tmp.q <- list(tmp.q)
            Rwk <- c(Rwk, tmp.q)
            tmp.q <- rkEval(ssr.obj$expand.call$rk[[i]], ssr.obj$data[[i]], 
                dataInit[[i]])
            if (!is.list(tmp.q)) 
                tmp.q <- list(tmp.q)
            thetaUsed <- thetaUsed[-c(1:length(tmp.q))]
            rwk <- c(rwk, lapply(tmp.q, function(x, z1) z1 * 
                x, z1 = ssr.obj$delta1[, i]))
        }
    
    phi.new <- NULL
    if (!is.null(terms1)) {
        for (i in 1:ncol(swk)) {
            phi.new <- cbind(phi.new, as.vector(phi[, i] %*% 
                t(as.matrix(terms1[, i]))))
        }
    }
    if (is.matrix(rwk)) 
        rwk <- list(rwk)
    if (is.matrix(Rwk)) 
        Rwk <- list(Rwk)
    qSum <- Rsum <- 0
    xi <- xi2 <- 0
    for (i in 1:length(theta)) {
        rwk[[i]] <- rwk[[i]] * theta[i]
        xi2 <- xi2 + as.vector(rwk[[i]]) %*% t(as.matrix(terms2[, 
            i]))
        Rsum <- Rsum + (as.vector(Rwk[[i]]) * theta[i]) %*% t(as.matrix(terms2[, 
            i]))
        qSum <- qSum + qwk[[i]] * theta[i]
        if (!is.null(ssr.obj$weight)) 
            rwk[[i]] <- ssr.obj$weight %*% rwk[[i]]
        xi <- xi + as.vector(rwk[[i]]) %*% t(as.matrix(terms2[, 
            i]))
    }
    y <- ssr.obj$y
#    if (is.dmudr && (ssr.obj$control$algorithm == "a2")) {
    if(is.dmudr){
        dsidr.fit <- dsidr(s = swk, q = qSum, y = y, weight = ssr.obj$weight, 
            vmu = ssr.obj$rkpk.obj$vmu, tol = ssr.obj$control$rkpk.control$tol, 
            job = ssr.obj$control$rkpk.control$job, limnla = ssr.obj$rkpk.control$limnla)
    }
    else dsidr.fit <- ssr.obj$rkpk.obj
    if (pstd) {
        b <- dsidr.fit$varht/(10^dsidr.fit$nlaht)
        if (nnull > 0) 
            dsmsV <- matrix(dsms(dsidr.fit)$sms, ncol = nnull)
        else dsmsV <- 0
        cr <- 0
        dr <- 0
        for (i in 1:length(rwk)) {
            dcrdrV <- dcrdr(dsidr.fit, rwk[[i]])
            cr <- cr + dcrdrV$cr %*% t(as.matrix(terms2[, i]))
            if (nnull > 0) 
                dr <- dr + dcrdrV$dr %*% t(as.matrix(terms2[, 
                  i]))
        }
        if (nnull > 0) {
            ciMain <- apply(phi, 1, function(x, terms1, w) diag(terms1 %*% 
                (x * w) %*% (x * t(terms1))), terms1 = terms1, 
                w = dsmsV)
            ciMain <- matrix(as.vector(ciMain), nrow = nrow(terms), 
                byrow = FALSE)
        }
        ciRK <- matrix(as.vector(cr) * as.vector(xi), ncol = nobs, 
            byrow = TRUE)
        ciRK <- matrix(apply(ciRK, 1, sum), nrow = eval.len, 
            byrow = FALSE)
        Rsum <- apply(Rsum, 2, function(x, r) x[seq(1, length(x), 
            length = r)], eval.len)
        ciRK <- Rsum - ciRK
        if (nnull > 0) {
            ciCross <- matrix(as.vector(t(phi.new)) * as.vector(dr), 
                ncol = nnull, byrow = TRUE)
            ciCross <- matrix(apply(ciCross, 1, sum), nrow = eval.len, 
                byrow = FALSE)
            ci <- t(ciMain) + ciRK - 2 * ciCross
        }
        else ci <- ciRK
        ci <- ci * (ci > 0)
    }
    ccc <- object$funcCoef$c
    ddd<- object$funcCoef$d
    fit <- apply(xi2, 2, function(x, w, r) matVec.prod(matrix(x, 
        ncol = r, byrow = TRUE), as.vector(w), left = FALSE), 
        w = ccc, r = nobs)
    if (nnull > 0) 
        fit <- fit + phi %*% matDiag.prod(t(terms1), ddd)
    if (pstd) {
        if (ncol(ci) == 1) 
            ci <- as.vector(ci)
        resul <- list(fit = fit, pstd = sqrt(ci * b))
    }
    else resul <- list(fit = fit, pstd = NULL)
    class(resul) <- c("bCI", "list")
    resul
}

## end of update

snm.control<- function(rkpk.control=list(job = -1, tol = 0, init = 0, limnla=c(-10, 3),  
              varht=NULL, theta= NULL, prec = 1e-006, maxit = 30),   
              nlme.control= list(returnObject=TRUE, maxIter=5), 
              prec.out=0.0005, maxit.out=30, converg= "COEF", incDelta=0.001) 
{ 
 
 ctrlvals.rk<-list(job = -1, tol = 0, init = 0, limnla=c(-10, 3),  
              varht=NULL, theta= NULL, prec = 1e-006, maxit = 30) 
 if(!(missing(rkpk.control) || is.null(rkpk.control))){ 
      for(nam in names(rkpk.control)) ctrlvals.rk[[nam]]<-rkpk.control[[nam]] 
 } 
 ctrlvals.nlme<- list(returnObject=TRUE, maxIter=5)        
 
 if(!(missing(nlme.control) || is.null(nlme.control))){ 
      for(nam in names(nlme.control)) ctrlvals.nlme[[nam]]<-nlme.control[[nam]] 
 } 
 if(missing(prec.out) || is.null(prec.out)) prec.out<- 0.0005 
 if(missing(maxit.out) || is.null(maxit.out)) maxit.out<- 30 
 if(missing(converg) || is.null(converg)) converg<- "COEF" 
 if(missing(incDelta) || is.null(incDelta)) incDelta<- 0.001 
 
 list(rkpk.control=ctrlvals.rk, nlme.control=ctrlvals.nlme, prec.out=prec.out, 
     maxit.out=maxit.out, converg=converg, incDelta=incDelta) 
} 

print.snm<-
function(x, ...)
{
   dd <- x$nlmeObj$dims
   cat("Semi-parametric nonlinear mixed-effects model fit by ")
   cat(ifelse(x$nlmeObj$method == "REML", "REML\n", "maximum likelihood\n"))

   cat("  Model:", deparse(as.vector(x$call$formula)), "\n")
   cat("  Data:", deparse(x$call$data), "\n")
   cat("  Log-", ifelse(x$nlmeObj$method == "REML", "restricted-", ""), 
		"likelihood: ", format(x$nlmeObj$logLik), "\n", sep = "")   
   cat("\n")
  fixF <- x$call$fixed
  if(inherits(fixF, "formula") || is.call(fixF)) {
		cat("Fixed:", deparse(as.vector(x$call$fixed)), "\n")
	}
	else {
		cat("Fixed:", deparse(lapply(fixF, function(el)
		as.name(deparse(as.vector(el))))), "\n")
	}
  print(fixef(x$nlmeObj))
  cat("\n")

  print(summary(x$nlmeObj$modelStruct), sigma=x$nlmeObj$sigma)
  
  switch(x$vmu,
          "v"= cat("GCV"),
          "m"= cat("GML"),
          "u"= cat("UBR"))
  cat(" estimate(s) of smoothing parameter(s): ")
  if(is.null(x$theta)) cat(format(10^x$lambda), "\n")
  else cat(format(10^x$theta), "\n")
  
  cat("Equivalent Degrees of Freedom (DF): ", format(x$funcDf), "\n")

  cat("\nNumber of Observations:", dd[["N"]]) 
  cat("\nNumber of Groups: ")
  Ngrps <- dd$ngrps[1:dd$Q]
  if((lNgrps <- length(Ngrps)) == 1) {
# single nesting
		cat(Ngrps, "\n")
	}
	else {
# multiple nesting
		sNgrps <- 1:lNgrps
		aux <- rep(names(Ngrps), sNgrps)
		aux <- split(aux, array(rep(sNgrps, lNgrps), c(lNgrps, lNgrps))[
			!lower.tri(diag(lNgrps))])
		names(Ngrps) <- unlist(lapply(aux, paste, collapse = " %in% "))
		cat("\n")
		print(rev(Ngrps))
	}
}
summary.snm<-function(object, ...){
## summary ob snm fit
  resul<- list(call=object$call, 
                vmu=object$vmu,
                nlme.fit=summary(object$nlmeObj),
                times=object$forCI$control$maxit.out,
                iter=object$iteration)
  resul$nlme.fit$call$fixed<-object$call$fixed
  resul$nlme.fit$call$random<-object$call$random
  resul$nlme.fit$call$model<- object$call$formula
  if(is.null(object$theta)) resul$lamb<- 10^object$lambda
  else resul$lamb<- 1/(10^object$theta)    
  resul$df<- object$funcDf
  class(resul)<- "summary.snm" 
  resul
 }
print.summary.snm<- function(x, ...){
   dd <- x$nlme.fit$dims
   cat("Semi-parametric Nonlinear Mixed Effects Model fit\n")
   cat("  Model:", deparse(as.vector(x$call$formula)), "\n")
   cat("  Data:", deparse(x$call$data), "\n")
   cat("\n")
## summary of nlme.fit 
   print(data.frame(AIC = x$nlme.fit$AIC+2*x$df, 
          BIC = x$nlme.fit$BIC+x$df*log(x$nlme.fit$dims$N), 
          logLik = x$nlme.fit$logLik, row.names = " "))  

  cat("\n")
  print(summary(x$nlme.fit$modelStruct), sigma = x$nlme.fit$sigma, 
        reEstimates = x$nlme.fit$coef$random)

  cat("Fixed effects: ")
	fixF <- x$nlme.fit$call$fixed
	if(inherits(fixF, "formula") || is.call(fixF)) {
		cat(deparse(as.vector(x$nlme.fit$call$fixed)), "\n")
	}
	else {
		cat(deparse(lapply(fixF, function(el)
		as.name(deparse(as.vector(el))))), "\n")
	}

## fixed effects t-table and correlations
  xtTab <- as.data.frame(x$nlme.fit$tTable)
  wchPval <- match("p-value", names(xtTab))
  for(i in names(xtTab)[ - wchPval]) {
     xtTab[, i] <- format(zapsmall(xtTab[, i]))
     }
  xtTab[, wchPval] <- format(round(xtTab[, wchPval], 4))
  if(any(wchLv <- (as.double(levels(xtTab[, wchPval])) == 0))) {
		levels(xtTab[, wchPval])[wchLv] <- "<.0001"
	}
  row.names(xtTab) <- dimnames(x$nlme.fit$tTable)[[1]]
  print(xtTab)
  if(nrow(x$nlme.fit$tTable) > 1) {
	corr <- x$nlme.fit$corFixed
	class(corr) <- "correlation"
	print(corr, title = " Correlation:")
   }

## summary of f
    switch(x$vmu,
          "v"= cat("\nGCV"),
          "m"= cat("\nGML"),
          "u"= cat("\nUBR"))

  cat(" estimate(s) of smoothing parameter(s):", format(x$lamb), "\n")  
  cat("Equivalent Degrees of Freedom (DF): ", format(x$df), "\n")
  if(x$times> x$iter[[1]])
       cat("\nConverged after", x$iter[[1]], "iterations\n")
  else cat("\nNot converged after", x$iter[[1]], "iterations\n")
}

predict.snm<- function(object, newdata=NULL, ...)
{
## To calculate prediction for snm object
  if(!inherits(object, "snm")) 
       stop("Input doesn't match")
     
  ssr.obj<- object$forCI  
 
  if(is.null(theta<- object$theta)) theta<- 1
  else theta<- 10^theta    

  nobs<- nrow(ssr.obj$s)
  nnull<- ssr.obj$rkpk.obj$nnull
  if(is.null(newdata)) eval.len<- nobs
  else   eval.len<- nrow(newdata)

## construct f's
  f.info<-getFunInfo(eval(object$call$func))
  funcName<- f.info$fName
  funcArg<- f.info$f.arg

  lengthF<- length(funcName)
#  assign("funcName", funcName, frame=1)
  spline.f0<- list(length(funcName))
  for(numF in 1:lengthF){            
    spline.f0[[numF]]<-function(...) {
      thisFunName<- as.character(match.call()[[1]])
      thisFunOrder<-match(thisFunName, funcName)
      res<-cbind(...)
      res.n<- nrow(res)
      res<- cbind(matrix(rep(0, (1+2*(thisFunOrder-1))*res.n), nrow=res.n), rep(1, res.n), res)
      if(length(funcName)== thisFunOrder) res 
      else
      cbind(res, matrix(rep(0, res.n*2*(length(funcName)-thisFunOrder)), nrow=res.n)) 
     }
  }   

  if(!missing(newdata)){
## order newdata
   grp<-getGroups(newdata, eval(object$call$groups))
   if(!is.null(object$call$data)){
     grpOld<- getGroups(eval(object$call$data), eval(object$call$groups))
    }
   else{
       grpName<- matrix(all.vars(eval(object$call$groups)[[2]]), nrow=1)
       grpOld<- data.frame(apply(grpName, 2, get))
       names(grpOld)<- grpName
   }
   if(inherits(grp, "factor")){
#      ordered(grp)<- levels(grpOld)
      grp<-ordered(grp,level=levels(grpOld))
      ord.ind<-order(grp)
      grp<- sort(grp)
      nobs.in <- list(table(grp))
      nobs.gp <- length(nobs.in[[1]])
     }
  else{
     for(numCol in 1:ncol(grpOld)){
#        ordered(grp[,numCol])<- levels(grpOld[,numCol])
         grp[,numCol]<-ordered(grp[,numCol], level=levels(grpOld[,numCol]))
       }        
       
      ord.ind<- do.call("order", grp) 
#      grp<-apply(grp, 2, order)
      nobs.in<- sapply(grp, function(x) table(sort(x)))
#      nobs.in<- rev(nobs.in)
      nobs.gp<- unlist(lapply(nobs.in, length))
   }    
   
     ord.org<- as.character(1:nrow(newdata))
     ord.org<- match(ord.org, ord.org[ord.ind])

## extract effect names
## fixed effects
  fixModelInfo<- getParaModelMatrix(eval(object$call$fixed), newdata, nrow(newdata))
  matrix.fix<- fixModelInfo$matrix.para
  length.fix<- fixModelInfo$length.para
  nameFixed<- fixModelInfo$para.name
  matrix.fix<- matrix.fix[ord.ind, , drop=FALSE]

##random effects  
  random<- reStruct(eval(object$call$random))
  re.form<- (formula(random))[[1]]
  nameRandom<- unlist(lapply(re.form, function(x) as.character(getResponseFormula(x)[[2]])))
  dataTmp<- newdata
  dataTmp[nameRandom]<-rep(1, nrow(dataTmp))
  matrix.random<- model.matrix(random, dataTmp)

  length.random<- unlist(lapply(re.form, 
          function(x, z) ncol(model.matrix2(getCovariateFormula(x), z)), z=dataTmp))
  matrix.random<- matrix.random[ord.ind, , drop=FALSE]

## combine fixed and random effects
  paraName<- c(nameFixed, nameRandom[match(nameRandom, nameFixed, nomatch=0)==0])

  coef.nlme<- object$nlmeObj$coef
  start.fixed <- as.vector(coef.nlme$fixed)
  start.random<-coef.nlme$random

## prepare for the evaluations of 3 delta's 
 
   para.v<- getParaValue(length.fix, length.random, matrix.fix, 
        matrix.random,start.fixed, start.random, paraName, nameRandom, 
         nrow(newdata), nameFixed, nobs.in)

    dataAug<- data.frame(newdata, para.v[ord.org,, drop=FALSE])

     for(i in 1:lengthF){
        assign(funcName[i], spline.f0[[i]])
     }

    delta<- eval(ssr.obj$expand.call$formula[[3]], dataAug)
    delta0<- delta[,1]
    delta<- apply(delta[,-1], 2, function(x, x0) x-x0, x0=delta0)
    dataInit<- list(lengthF)
    delta1<- NULL 
    fArgLen<- c(0, cumsum(unlist(lapply(funcArg, length))+1))
    for(i in 1:lengthF){
       dataTmp<- delta[, (fArgLen[i]+1):(fArgLen[i+1])]
       delta1<- cbind(delta1, dataTmp[,1])
       dataTmp<- dataTmp[,-1]/dataTmp[,1]
       dataInit[[i]]<- data.frame(dataTmp)
       names(dataInit[[i]])<- funcArg[[i]]
      }

     fit.y<- delta0

     oldCoefD<- object$funcCoef$d
     thetaUsed<- theta
     for(i in 1:lengthF){
        if(length(ssr.obj$expand.call$formF)>0 && !is.null(ssr.obj$expand.call$formF[[i]])){
          tmp.s<- model.matrix2(ssr.obj$expand.call$formF[[i]], dataInit[[i]]) 
          fit.y<- fit.y+as.vector((delta1[,i]*tmp.s)%*%oldCoefD[1:ncol(tmp.s)])
          oldCoefD<- oldCoefD[-c(1:ncol(tmp.s))]
        }
      
        tmp.q<- rkEval(ssr.obj$expand.call$rk[[i]], ssr.obj$data[[i]], dataInit[[i]])
        if(!is.list(tmp.q)) tmp.q<- list(tmp.q)
        fit.y<- fit.y+ matVec.prod(sumList(list.prod(tmp.q,  thetaUsed[1:length(tmp.q)])), 
                   object$funcCoef$c*ssr.obj$delta1[, i])*delta1[,i]
        thetaUsed<- thetaUsed[-c(1:length(tmp.q))]
     }
   }
# calculate fits
  else fit.y<- as.vector(object$fitted)
  fit.y
}    

plot.bCI<-
function(x, x.val = NULL, type.name = NULL, ...)
{
        fit<- x
	if(match("listbCI", class(fit), nomatch=FALSE)){
              f.name<- names(fit)
              for(obj in 1:length(fit)) print(plot.bCI(fit[[obj]], x.val=x.val, type.name=f.name[obj], ...))
        }
	else{
        num <- dim(fit[[1]])

        if(is.null(num))
                num <- c(length(fit[[1]]), 1)
        if(!is.null(x.val)) {
                if(length(x.val) != num[1])
                        stop("length of x.val not match")
        }
        else x.val <- c(1:num[1])
        pt.est <- as.vector(fit$fit)
        pstd <- fit$pstd
        if(is.null(type.name))
                type.name <- paste("fit", 1:num[2], sep = "")
	else if(length(type.name)<num[2]) type.name<- paste(type.name, 1:num[2], sep="")
        pt.dat <- data.frame(estimate = as.vector(fit$fit), x = rep(x.val, num[
                2]), f.type = rep(type.name, rep(num[1], num[2])))
        if(!is.null(pstd)) {
                up.dat <- data.frame(estimate = as.vector(fit$fit) + 1.96 * 
                        as.vector(pstd), x = rep(x.val, num[2]), f.type = rep(
                        type.name, rep(num[1], num[2])))
                low.dat <- data.frame(estimate = as.vector(fit$fit) - 1.96 * 
                        as.vector(pstd), x = rep(x.val, num[2]), f.type = rep(
                        type.name, rep(num[1], num[2])))
        }
        if(is.null(pstd)) {
                resul <- xyplot(estimate ~ x | f.type, data = pt.dat, ...)
        }
        else {
                resul <- xyplot2(estimate ~ x | f.type, data = list(pt.dat, 
                        up.dat, low.dat), scale.v = "free", ...)
        }
        resul
}}

plot.ssr<-
function(x, ask = FALSE, ...)
{
   ssr.obj<- x
## plot on ssr object
	choices <- c("All the following plots", "Estimate of function with CIs",
		"Residuals against Fitted values", 
		"Response against Fitted values", "Normal QQplot of Residuals")
	choices <- substring(choices, 1, 40)
	tmenu <- paste("", choices)
	pick <- 2
	ask.now <- ask
	while(pick <= length(tmenu) + 2) {
		if(ask.now)
			pick <- menu(tmenu, title = 
				"\nMake a plot selection (or 0 to exit):\n") + 
				1
		switch(pick,
			invisible(return(ssr.obj)),                       
			{
				ask.now <- FALSE
			}
			,
## fitted values
			{
## calculate fits and STDs
				ci <- predict(ssr.obj)
				if(ssr.obj$family != "gaussian") {
				  switch(ssr.obj$family,
				    binary = ci <- lapply(ci, alogit),
				    binomial = {
				      ci <- lapply(ci, alogit)
				    }
				    ,
				    poisson = ci <- lapply(ci, exp),
				    gamma = ci <- lapply(ci, exp))
				}
				if((length(ssr.obj$q) == 1) && (length(ssr.obj$
				  call$rk[[2]]) == 1)) {
				  t.data <- ssr.obj$data
				  xarg <- as.character(ssr.obj$call$rk[[2]])
				  if(is.data.frame(t.data))
				    xarg <- t.data[[xarg]]
#				  else xarg <- eval(parse(text = xarg), local = 0)
                                  else xarg<- eval(parse(text=xarg))
				}
				tempx <- 1:length(ci$fit)
				if((length(ssr.obj$q) == 1) && (length(ssr.obj$
				  call$rk[[2]]) == 1)) {
				  tempx <- sort(xarg)
				  ci[[1]] <- ci[[1]][order(xarg)]
				  ci[[2]] <- ci[[2]][order(xarg)]
				}
				up <- as.vector(ci$fit) + 1.96 * as.vector(ci$
				  pstd)
				down <- as.vector(ci$fit) - 1.96 * as.vector(ci$
				  pstd)
				if((length(ssr.obj$q) > 1) || (length(ssr.obj$
				  call$rk[[2]]) > 1))
				  x.char <- "x"
				else x.char <- as.character(ssr.obj$call$rk[[2
				    ]])
				plot(tempx, as.vector(ci$fit), xlab = x.char, 
				  ylab = "Spline Estimate of f", ylim = c(min(
				  down), max(up)), type = "l")
				lines(tempx, up, lty = 2)
				lines(tempx, down, lty = 2)
			}
			,
##resi vs fitted
			{
				plot(ssr.obj$fit, ssr.obj$resi, xlab = 
				  "Fitted Values", ylab = "Residuals")
			}
			,
##response vs fitted
			{
				if(ssr.obj$family != "binomial")
				  y.tmp <- ssr.obj$y
				else y.tmp <- ssr.obj$y[, 2]/ssr.obj$y[, 1]
				plot(ssr.obj$fit, y.tmp, xlab = "Fitted values",
				  ylab = "Response")
			}
			,
## qqplot
			{
				qqnorm(ssr.obj$residuals, ylab = "Residuals", 
				  xlab = "Quantiles of Standard Normal")
				qqline(ssr.obj$residuals)
			}
			)
		if(!ask.now)
			pick <- pick + 1
		if(pick == length(tmenu) + 2)
			ask.now <- ask
	}
	invisible()
}

### some addtional functions
"count.rows"<-
function(x, ...)
UseMethod("count.rows")

"count.rows.default"<-
function(x, ...)
{
	if(is.matrix(x))
		count.rows.matrix(x, ...)
	else if(is.list(x))
		max(sapply(x, length))
	else length(x)
}

"select.col"<-
function(x, column)
UseMethod("select.col")

"select.col.default"<-
function(x, column)
{
	if(is.matrix(x))
		x[, column]
	else if(is.list(x))
		x[[column]]
	else x
}

"order.rows"<-
function(target, columns.to.sort.by)
{
	ret <- seq(length = count.rows(target))
	runs <- list(ret)
	run.count <- ifelse(length(ret) > 0, 1, 0)
	run.index <- 1
	by.count <- length(columns.to.sort.by)
	by.index <- 1
	by.last.run <- run.count
	while(run.index <= run.count) {
		this.run <- runs[[run.index]]
		this.ret <- ret[this.run]
		values.to.sort.by <- select.col(target, columns.to.sort.by[
			by.index])[this.ret]
		values.order <- order(values.to.sort.by)
		ret[this.run] <- this.ret[values.order]
		if(by.index < by.count) {
			run.sum <- 0
			for(run.check in rle(values.to.sort.by[values.order])$
				lengths) {
				run.sum <- run.sum + run.check
				if(run.check > 1) {
				  run.count <- run.count + 1
				  runs[[run.count]] <- this.run[seq(to = 
				    run.sum, length = run.check)]
				}
			}
			if(run.index == by.last.run) {
				by.index <- by.index + 1
				by.last.run <- run.count
			}
		}
		run.index <- run.index + 1
	}
	ret
}

"count.rows.matrix"<-
function(x, ...)
dim(x)[1]

inc <- function(y, x, spar="v", limnla=c(-6,0), grid=x, 
                prec=1.e-6, maxit=50, verbose=F)
{
 n <- length(x)
 org.ord <- match(1:n, (1:n)[order(x)])
 s.x <- sort(x)
 s.y <- y[order(x)]
 x1 <- c(0, s.x[-n])
 x2 <- s.x-x1
 q.x <- as.vector(rep(1,3)%o%x1+c(0.1127017,0.5,0.8872983)%o%x2)

 # function for computing derivatives
 k1 <- function(x) x-.5
 k2 <- function(x) ((x-.5)^2-1/12)/2
 dk2 <- function(x) x-.5
 dk4 <- function(x) sign(x)*((abs(x)-.5)^3/6-(abs(x)-.5)/24)
 drkcub <- function(x,z) dk2(x)%o%k2(z)-dk4(x%o%rep(1,length(z))-rep(1,length(x))%o%z)

 # compute starting value
 ini.fit <- ssr(s.y~I(s.x-.5), cubic(s.x))
 g.der <- ini.fit$coef$d[2]+drkcub(q.x,x)%*%ini.fit$coef$c
 h.new <- abs(g.der)+0.005

 # begin iteration
 iter <- cover <- 1  
 h.old <- h.new 
 repeat {
  if(verbose) cat("\n Iteration: ", iter)
  yhat <- s.y-int.s(h.new*(1-log(h.new)),s.x)        
  smat <- cbind(int.s(h.new, s.x), int.s(h.new*q.x,s.x))
  qmat <- int1(s.x,h.new)
  fit <- ssr(yhat~smat, qmat, spar=spar, limnla=limnla)
  if(verbose) cat("\nSmoothing parameter: ", fit$rkpk$nlaht)
  dd <- fit$coef$d
  cc <- as.vector(fit$coef$c)
  h.new <- as.vector(exp(cc%*%int.f(s.x,q.x,h.new)+dd[2]+dd[3]*q.x))
  cover <- mean((h.new-h.old)^2)     
  h.old <- h.new
  if(verbose) cat("\nConvergent Criterion: ", cover, "\n")
  if(cover<prec || iter>(maxit-1)) break
  iter <- iter + 1 
 }
 if(iter>=maxit) print("convergence not achieved!")
 y.fit <- (smat[,1]+fit$rkpk$d[1])[org.ord]
 f.fit <- as.vector(cc%*%int.f(s.x,grid,h.new)+dd[2]+dd[3]*grid)
 x1 <- c(0, grid[-length(grid)])
 x2 <- grid-x1
 q.x <- as.vector(rep(1,3)%o%x1+c(0.1127017,0.5,0.8872983)%o%x2)
 h.new <- as.vector(exp(cc%*%int.f(s.x,q.x,h.new)+dd[2]+dd[3]*q.x))
 y.pre <- int.s(h.new,grid)+fit$rkpk$d[1]
 sigma <- sqrt(sum((y-y.fit)^2)/(length(y)-fit$df))
# f.fit <- as.vector(cc%*%int.f(s.x,grid,h.new)+dd[2]+dd[3]*grid)
# f.pstd <- predict(fit,terms=c(0,1,1,1))$pstd
 list(fit=fit, iter=c(iter, cover), pred=list(x=grid,y=y.pre,f=f.fit), 
      y.fit=y.fit, sigma=sigma)
}    

# functions for computing integrals
int.s <- function(f, x, low = 0) {
 n <- length(x); x <- c(low, x)
 .C("integral_s", as.double(f), as.double(x),
  as.integer(n), val = double(n))$val
}

int.f <- function(x, y, f, low = 0) {
 nx <- length(x); ny <- length(y); x <- c(low, x)
 res <- .C("integral_f", as.double(x),
  as.double(y), as.double(f), as.integer(nx),
  as.integer(ny), val=double(nx*ny))$val
 matrix(res, ncol=ny, byrow=F)
}

int1 <- function(x, f.val, low = 0) {
 n <- length(x); x <- c(low, x)
 if(length(f.val) != 3 * n) stop("input not match")
 res <- matrix(.C("integral_1", as.double(x),
  as.double(x), as.double(f.val), as.integer(n),
  as.integer(n), val = double(n * n))$val, ncol = n)
 apply(res, 1, cumsum)
}


#.onLoad<- function(...){
  #require(nlme)
  #library.dynam("assist")
#}
