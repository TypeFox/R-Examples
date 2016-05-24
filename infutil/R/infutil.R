require(ltm)

infoGPCM <- ltm:::infoGPCM
probs <- ltm:::probs

#################################################
### Functions for calculating utility, regret ###
#################################################

## count.rows function from Claudia Beleites; http://rwiki.sciviews.org/doku.php?id=tips:data-frames:count_and_extract_unique_rows ##

count.rows <- function(x) { 
	if (is.matrix (x) && (dim (x) [2] == 1))
	x <- as.vector (x) 
 
	order.x <- do.call(order,as.data.frame(x))
   
	if (is.vector (x)) {
		equal.to.previous <-
		x[tail(order.x,-1)] == x[head(order.x,-1)]
		} else {
			equal.to.previous <-
			rowSums(x[tail(order.x,-1),] != x[head(order.x,-1),])==0
			}	
 
	indices <-  split (order.x, cumsum (c (TRUE, !equal.to.previous)))
 
	if (is.vector (x)) {
	x <- x [sapply (indices, function (x) x [[1]]), drop = FALSE]
		} else {
			x <- x [sapply (indices, function (x) x [[1]]), ,drop = FALSE]
			}
  
	data.frame (counts = sapply (indices, length) , ind = I(indices), x) 
	}

	
## infoprobs, cprobs, and iprobs originally by Dimitris Rizopoulos for the ltm package ##
## they are modified here to accept scalar latent trait values (z)  ##
	
infoprobs <- function (betas, z)	{
	cpr <- cprobs(betas, z)
	ipr <- iprobs(betas, z)
	sum.cprs <- lapply(cpr, function(x) {
	    nr <- nrow(x)
	    if(ncol(x) == 1){
	    	t((1 - c(x[1, ], x[-nr, ] + x[-1, ]))^2)
	    	} else if(ncol(x) > 1){
	    		t((1 - rbind(x[1, ], x[-nr, ] + x[-1, ]))^2)
	    		}
	    		
	    })

	betas. <- sapply(betas, function(x) x[length(x)])

	for (i in 1:length(betas)) sum.cprs[[i]] <- betas.[i]^2 * ipr[[i]] * sum.cprs[[i]]

	do.call(cbind, lapply(sum.cprs, rowSums))
	
	}	
	
cprobs <- function (betas, z, eps = .Machine$double.eps^(1/3)){
    lapply(betas, function(x, z) {
        nx <- length(x)
        out <- plogis(x[-nx] - matrix(x[nx] * z, nx - 1, length(z), 
            TRUE))
        if (any(ind <- out == 1)) 
            out[ind] <- 1 - eps
        if (any(ind <- out == 0)) 
            out[ind] <- eps
        rbind(out, 1)
    }, z = z)
    }

iprobs =
function (betas, z) 
{
    n <- length(z)
    gammas <- lapply(betas, function(x) {
        nx <- length(x)
        cbind(plogis(matrix(x[-nx], n, nx - 1, TRUE) - x[nx] * z), 1)
    })
    lapply(gammas, function(x) {
        nc <- ncol(x)
        if(n==1){
          c(x[, 1], x[, 2:nc] - x[, 1:(nc - 1)])
          } else if(n>1){
          	cbind(x[, 1], x[, 2:nc] - x[, 1:(nc - 1)])
          	}
    })
}	



## inf.func.ltm and logL.ltm contain code originally by Dimitris Rizopoulos for the ltm package ##

inf.func.ltm <- 
function (object, items = NULL, ...) 
{
    if (!class(object) %in% c("grm", "gpcm", "ltm", "rasch", 
        "tpm")) 
        stop("'object' must inherit from either class 'grm', class 'gpcm', class 'ltm', class 'rasch' or class 'tpm'.\n")
    p <- ncol(object$X)
    itms <- if (!is.null(items)) {
        if (!is.numeric(items) && length(items) > p) 
            stop("'items' should be a numeric vector of maximum length ", 
                p, ".\n")
        if (any(!items %in% 1:p)) 
            stop("'items' should contain numbers in: ", paste(1:p, 
                collapse = ", "), " indicating the items.\n")
        items
    }
    else 1:p
    if (class(object) == "ltm" && (object$ltst$factors > 1 | 
        any(unlist(object$ltst[c("inter", "quad.z1", "quad.z2")])))) 
        stop("Information is currently computed only for the one-factor two-parameter logistic model.\n")
    f <- function(z) {
        switch(class(object),
            grm =  rowSums(infoprobs(object$coefficients, z)[, itms, drop = FALSE]),
            gpcm = rowSums(infoGPCM(object$coefficients, z, object$IRT.param)[, itms, drop = FALSE]),
            ltm = {
            betas <- object$coefficients
            Z <- cbind(1, z)
            mat <- t(t(plogis(Z %*% t(betas)) * (1 - plogis(Z %*% t(betas)))) * betas[, 2]^2)
            rowSums(mat[, itms, drop = FALSE])
        }, rasch = {
            betas <- object$coefficients
            Z <- cbind(1, z)
            mat <- betas[1, 2]^2 * plogis(Z %*% t(betas)) * (1 - plogis(Z %*% t(betas)))
            rowSums(mat[, itms, drop = FALSE])
        }, tpm = {
            thetas <- object$coefficients
            Z <- cbind(1, z)
            betas <- thetas[, 2:3]
            cs <- plogis(thetas[, 1]) * object$max.guessing
            pi. <- plogis(Z %*% t(betas))
            cs <- matrix(cs, nrow(Z), p, TRUE)
            pi <- cs + (1 - cs) * pi.
            pqr <- pi * (1 - pi) * (pi./pi)^2
            mat <- t(t(pqr) * betas[, 2]^2)
            rowSums(mat[, itms, drop = FALSE])
        })
    }

    f
    
}



Jeffreys <- function(ltm.obj=NULL, inf.mat=NULL, inf.func=NULL, return="prior", spl.method="natural", range.int=c(-Inf, Inf)){

	if(sum(is.null(ltm.obj), is.null(inf.mat), is.null(inf.func)) != 2)
		stop("One of ltm.obj, inf.mat, or inf.func must be supplied.\n")
		
	int.lower <- range.int[1]
	int.upper <- range.int[2]

	if(((class(ltm.obj) == "grm") | (class(ltm.obj) == "gpcm")) & any(abs(range.int) > 50)){ # Inf values cause problems for grm and gpcm for some reason
		if(int.lower < -50) { int.lower <- -50 }
		if(int.upper >  50) { int.upper <-  50 }
		}
		
	if(!is.null(inf.func)){

		Jp.nc <- integrate(function(x){sqrt(inf.func(x))}, int.lower, int.upper)$val	
		
		} else if(!is.null(inf.mat)){
		
			inf.func.interp <- splinefun(inf.mat[,1], inf.mat[,2], method=spl.method)
			
			inf.func <- function(theta){
				out <- inf.func.interp(theta)
				out[out < 0] <- 0
				out
				}

			Jp.nc <- integrate(function(x){sqrt(inf.func(x))}, int.lower, int.upper)$val
				
			} else if(!is.null(ltm.obj)){

				inf.func <- inf.func.ltm(ltm.obj)

				Jp.nc <- integrate(function(x){sqrt(inf.func(x))}, int.lower, int.upper)$val
					
				}
		
	Jp <- function(theta){
		sqrt(inf.func(theta)) / Jp.nc
		}
		
	switch(return,
		prior = Jp,
		nc = Jp.nc,
		both = list(prior = Jp, nc = Jp.nc)	
		)
				
	}

	
rJeffreys <- function(n, prior, range.int=c(-Inf, Inf)){

	int.lower = range.int[1]
	int.upper = range.int[2]

	Jpinv.min <- function(xval, prob){ abs(prob - integrate(prior, int.lower, xval)$val) }

	Jpinv <- function(unvar){ nlminb(qnorm(unvar), Jpinv.min, prob=unvar, lower=int.lower, upper=int.upper)$par }

	sapply(runif(n), Jpinv)
	
	}
	
	

logL.ltm <- function(z, dat, i, ltm.obj){

	switch(class(ltm.obj),
	
		ltm = { logL <- function (z, dat, ltm.obj, i) {
				betas = ltm.obj$coef		
				y = dat[i, ]	
				Z <- c(1, z)
				names(Z) <- c("(Intercept)", "z1")
				Z <- Z[match(colnames(betas), names(Z))]
				pr <- probs(c(betas %*% Z))
				sum(dbinom(y, 1, pr, log = TRUE), na.rm = TRUE)
				}
		      },
		      
		rasch = { logL <- function (z, dat, ltm.obj, i) {
				betas = ltm.obj$coef
				y = dat[i, ]
				pr <- probs(c(betas %*% c(1, z)))
				sum(dbinom(y, 1, pr, log = TRUE), na.rm = TRUE)
				}
		        },
		        
		tpm = { logL <- function (z, dat, ltm.obj, i) {
				thetas = ltm.obj$coef
				betas <- thetas[, 2:3]
				y = dat[i, ]			
				cs <- plogis(thetas[, 1]) * ltm.obj$max.guessing
				pr <- cs + (1 - cs) * probs(c(betas %*% c(1, z)))
				sum(dbinom(y, 1, pr, log = TRUE), na.rm = TRUE)
				}
		      },
		      
		gpcm = { if(min(dat) < 1) {stop("Responses in data must be coded 1, 2, 3, ...")}
		
			 logL <- function (z, dat, ltm.obj, i) {
				betas = ltm.obj$coefficients
				p = length(betas)				
				y = dat[i, ]
				log.prs <- ltm:::crf.GPCM(betas, z, IRT.param = ltm.obj$IRT.param, log = TRUE)
				log.pxz <- numeric(p)
				
				for (j in 1:p) {
					log.pxz[j] <- if (!is.na(y[j])) log.prs[[j]][y[j]] else 0
					}
					
				if(any(!is.finite(log.pxz))){ 
					log.pxz[log.pxz == -Inf] = log(.Machine$double.xmin)
					log.pxz[log.pxz ==  Inf] = log(.Machine$double.xmax)
					} 
				
				sum(log.pxz, na.rm = TRUE)
				
				}      
		       },
		       
		grm = { if(min(dat) < 1) {stop("Responses in data must be coded 1, 2, 3, ...")}
		
			logL <- function (z, dat, ltm.obj, i) {
				betas = ltm.obj$coefficients
				p = length(betas)
				y = dat[i, ]
				gammas <- lapply(betas, function (x) {
					n <- length(z)
					nx <- length(x)
					c(plogis(matrix(x[-nx], n, nx - 1, TRUE) - x[nx] * z), 1)
					})
				log.prs <- lapply(gammas, function (x) {
					nc <- length(x)
					prs <- c(x[1], x[2:nc] - x[1:(nc - 1)])
					prs[prs == 0] = .Machine$double.eps
					log(prs)
					})
				log.pxz <- numeric(p)
				
				for (j in 1:p) {
					log.pxz[j] <- if (!is.na(y[j])) log.prs[[j]][y[j]] else 0
					}
					
				sum(log.pxz, na.rm = TRUE)
				}
		       }
		)	
	
	Vectorize(logL, "z")(z=z, dat=dat, i=i, ltm.obj=ltm.obj)
	
	
	
	}



iota <-
function(ltm.obj, logL.fun, fscore.obj=NULL, data=NULL, prior=NULL, theta0=NULL, range.int=c(-Inf, Inf), range.theta = c(-10, 10)){


	if( (!is.null(fscore.obj) & !is.null(data)) | (is.null(fscore.obj) & is.null(data)) ){ stop("One of fscore.obj or data must be supplied.\n") }
	if( !is.null(fscore.obj) & (class(fscore.obj) != "fscores") ){ stop("fscore.obj must be a fscores object.\n") }
	if( (class(prior) != "function") | (is.null(prior)) ){ stop("A prior must be supplied in the form of a function.\n") }
	if(any(!is.finite(range.theta))) {stop("range.theta must be finite.")}
		
	int.lower = range.int[1]
	int.upper = range.int[2]

	if(((class(ltm.obj) == "grm") | (class(ltm.obj) == "gpcm")) & any(abs(range.int) > 10)){ # Inf values cause problems for grm and gpcm for some reason
		if(int.lower < -10) { int.lower = -10 }
		if(int.upper >  10) { int.upper =  10 }
		}
	
	if( missing(logL.fun) & missing(ltm.obj) ) {stop("A log-likelihood function or ltm object must be supplied.\n")}	
	if( missing(logL.fun) & (!missing(ltm.obj)) ) {
	
		logL.fun <- function(z, dat, i){ logL.ltm(z=z, dat=dat, i=i, ltm.obj=ltm.obj) }
	
		}
	
	if(!is.null(data)){
		
		data.unique = count.rows(data)
		N = nrow(data)
		N.u = nrow(data.unique)
		p.xe = data.unique$counts / N
		dat.mat = as.matrix(data.unique[,-c(1,2)])
	
		} else if(!is.null(fscore.obj)){

			N = sum(fscore.obj$score.dat$Obs)
			N.u = nrow(fscore.obj$score.dat)
			p.xe = fscore.obj$score.dat$Obs / N
			dat.mat = as.matrix(fscore.obj$score.dat[,1:nrow(fscore.obj$coef)])
			
			}

	if( !is.null(theta0) ){
	
		if( (mode(theta0) == "numeric") & (length(theta0) > 1) & (length(theta0) != N.u) ) {
			stop("theta0 must be a scalar or equal in length to the number of unique response patterns.")
			} else if( mode(theta0) == "character" ) {
				if( length(theta0) > 1 ) {
					stop("theta0 must be a numeric scalar or vector, \"max.prior\", or NULL\n") 
					} else if( (theta0 != "max.prior") ) {
						stop("theta0 must be a numeric scalar or vector, \"max.prior\", or NULL\n") 
						}
				} 
	
		}
			
	
	
	
	I.x = vector(length=N.u)
	p.x = vector(length=N.u)

	dlogp.xt <- function(z, i){
			logL.fun(z=z, dat=dat.mat, i=i) 
			}

	dp.tx <- function(z, i){
			exp(dlogp.xt(z=z, i=i))*prior(z) 
			}
			
	dI.x <-  function(z, i, p.x.i){
			(dp.tx(z=z, i=i) / p.x.i) * (dlogp.xt(z=z, i=i) - log(p.x.i))			
			}
		
	for(i in 1:N.u)	{

		p.x[i] = integrate(dp.tx, int.lower, int.upper, i=i)$val
	        I.x[i] = integrate(dI.x, int.lower, int.upper, i=i, p.x.i=p.x[i])$val

		}

	I = sum(I.x*p.xe)
		
	if( (!is.null(theta0)) & (length(theta0) == 1) ){
		
		if(theta0 == "max.prior"){ theta0 = optimize(prior, range.theta, maximum=T)$max }
		
		logL0  = vector(length=N.u)
			
		for(i in 1:N.u)	{
        	
			logL0[i] = dlogp.xt(z=theta0, i=i)
					
			}
			
		} else if( (!is.null(theta0)) & (length(theta0) == N.u) ){
			
			logL0  = vector(length=N.u)
			
			for(i in 1:N.u)	{
        	
				logL0[i] = dlogp.xt(z=theta0[i], i=i)
					
				}
				
			}

	if( !is.null(data) ){
		
		I.x.exp = vector(length=N)
		p.x.exp = vector(length=N)
		
		for(i in 1:N.u){
			I.x.exp[data.unique$ind[[i]]] = I.x[i]		
			p.x.exp[data.unique$ind[[i]]] = p.x[i]
			}
		
		I.x = I.x.exp
		p.x = p.x.exp
			
		if( !is.null(theta0) ){
			logL0.exp = vector(length=N)
			for(i in 1:N.u){logL0.exp[data.unique$ind[[i]]] = logL0[i]}
			logL0 = logL0.exp
			}			
			
		}
		

	if( is.null(theta0) )  {
	
		list(I=I, I.x=I.x, p.x=p.x )
		
		} else if( !is.null(theta0) ) {
		
			list(I=I, I.x=I.x, p.x=p.x, logL0=logL0, logNL0.pval = 1-pchisq(-2*((logL0 - log(p.x))-I.x), df=1) )
			
			} 		
				
	}
	
	
	
iota.l <- function(x){
		if(class(x) %in% c("ltm", "grm", "gpcm")){
			nc = Jeffreys(x, return="nc")
			} else if( class(x) == "numeric" ) {
				nc = x
				} else { stop("x must be an ltm object or a normalizing constant.\n") }
			
		(1/2)*log(1/(2*pi*exp(1))) + log(nc)
			
		}

		
		
iota.u <- function(prior, range.int=c(-Inf, Inf)){
		if(class(prior) != "function"){ stop("prior must be a function.\n") }

		int.lower = range.int[1]
		int.upper = range.int[2]
		
		Hint <- function(theta){
			p.t = prior(theta)
			integrand = prior(theta)*log(prior(theta))
			
			integrand[p.t == 0] = 0
			
			integrand
			
			}
			
		-integrate(Hint, int.lower, int.upper)$val 
			
		}
	
		
		
nmru <- function(ltm.obj, range.int=c(-Inf, Inf)){

	int.lower = range.int[1]
	int.upper = range.int[2]

	if(!(class(ltm.obj) %in% c("ltm", "grm", "gpcm"))){ stop("x must be an ltm object.\n") }

	Jp.obj = Jeffreys(ltm.obj, return="both")
	
	prior = Jp.obj$prior
	nc    = Jp.obj$nc
		
	iota.lower = iota.l(nc)
	iota.upper = iota.u(prior, range.int=range.int)

	list(val=iota.lower/iota.upper, iota.l=iota.lower, iota.u=iota.upper)
	
	}
	
	
		
iota.c <-
function(ltm.obj, M=NULL, prior=NULL, logL.fun=NULL, rirm=NULL, range.int=c(-Inf, Inf)){

	if(is.null(M)){stop("The number of Monte Carlo replications M must be supplied.\n")}
	if((!is.null(prior)) & (class(prior) != "function") ){ stop("The reference prior must be supplied in the form of a function.\n") }
	
	int.lower = range.int[1]
	int.upper = range.int[2]
	
	if(((class(ltm.obj) == "grm") | (class(ltm.obj) == "gpcm")) & any(abs(range.int) > 10)){ # Inf values cause problems for grm and gpcm for some reason
		if(int.lower < -10) { int.lower = -10 }
		if(int.upper >  10) { int.upper =  10 }
		}

	if( missing(ltm.obj) & (is.null(prior) | is.null(logL.fun) | is.null(rirm)) ) {
		stop("A log-likelihood function, response generating function, and reference prior must be supplied in the absence of an ltm object.\n")
		}
		
	if( !missing(ltm.obj) ) {
		
		prior = Jeffreys(ltm.obj=ltm.obj)
		logL.fun <- function(z, dat, i){ logL.ltm(z=z, dat=dat, i=i, ltm.obj=ltm.obj) }
	
		switch(class(ltm.obj),
	        
		      ltm   = { rirm <- function(M, m.theta) { rmvlogis(n=M, thetas=coef(ltm.obj), IRT=ltm.obj$IRT.param, z.vals=m.theta) }  },
		      rasch = { rirm <- function(M, m.theta) { rmvlogis(n=M, thetas=coef(ltm.obj), IRT=ltm.obj$IRT.param, z.vals=m.theta) }  },
		      gpcm  = { rirm <- function(M, m.theta) { rmvordlogis(n=M, thetas=as.list(as.data.frame(t(coef(ltm.obj)))), IRT=ltm.obj$IRT.param, z.vals=m.theta, model=class(ltm.obj)) }  },
		      grm   = { rirm <- function(M, m.theta) { rmvordlogis(n=M, thetas=as.list(as.data.frame(t(coef(ltm.obj)))), IRT=ltm.obj$IRT.param, z.vals=m.theta, model=class(ltm.obj)) }  }

		      )
		      
		}
	
	m.theta = rJeffreys(M, prior, range.int=c(int.lower, int.upper))
		
	data = rirm(M=M, m.theta=m.theta)
		
	data.unique = count.rows(data)
	N.u = nrow(data.unique)
	dat.mat = as.matrix(data.unique[,-c(1,2)])

	
	p.x = vector(length=M)
	I.x = vector(length=M)
	
	dlogp.xt <- function(z, i){
			logL.fun(z=z, dat=dat.mat, i=i) 
			}

	dp.tx <- function(z, i){
			exp(dlogp.xt(z=z, i=i))*prior(z) 
			}
		
	for(i in 1:N.u)	{
        
		p.x[data.unique$ind[[i]]] = integrate(dp.tx, int.lower, int.upper, i=i)$val
		
		}
        
	for(i in 1:M){
		
		I.x[i] = logL.fun(z=m.theta[i], dat=data, i=i) - log(p.x[i])		
	
		}
		
	list(I = mean(I.x), se.I = sd(I.x) / sqrt(M) , I.x=I.x, p.x=p.x)
		
	}
		



