.packageName<-'sensory'

################
###Initialize###
################

combinewk <- function(weights=NULL, label=NULL)	{
	if (is.null(label)) stop('label is null')
	kw = label !=0
 	for (j in 1:ncol(weights)) weights[kw,j] = (label == j)[kw]
	return(weights)	
	}

ipar <- function(data,q=1) {
	p <- ncol(data)
	val <- list()
	val$mu <- rnorm(p, apply(data,2,mean,na.rm=TRUE),apply(data,2,sd, na.rm=TRUE)/nrow(data))
	S <- cov(data, use="pairwise.complete.obs")
	S[is.na(S)] = 0
	val$lambda <- matrix(eigen(S)$vectors[,1:q], nrow=p, ncol=q)
	val$tlambda <- t(val$lambda)
	val$psi <- diag(S - val$lambda %*% val$tlambda)
	val$psi[val$psi < 0] <- runif(sum(val$psi <0), 1/2, 1) 
	val$invS <- invSigma(tlam=val$tlambda, lam=val$lambda, psi=val$psi)
	return(val)
	}

rgpar <- function(data, g=2,q=1) {
	val <- list()
	for (k in 1:g) val[[k]] <- ipar(data=data,q=q)
	val$pi <- rdirichlet(n=1,rep(1,g))
	return(val)
	}	

igpar <- function(data, g=2, q=1) {
	val1 <- rgpar(data, g=1, q=q)
	lval <- as.numeric( kmeans( x=initial.ySg(data, par= val1[[1]], weights=NULL)$y, centers=g, iter.max=10)$cluster )
	valg <- rgpar(data, g=g, q=q)
	val <- EMn(data=data, igpar=valg, n=1, label=lval)$gpar
	return(val)
	}

#######################
###Woodbury Identity###
#######################

logdetSigma <- function(tlam=NULL,lam=NULL, psi=NULL, invS=NULL){
	lpsi <- sum(log(psi))
	if ( is.null(invS) ) invS <- invSigma(tlam=tlam,lam=lam, psi=psi)
	q <- ncol(lam)
	invT <- diag(rep(1,q)) - tlam %*% invS %*% lam
	loginvT <- log(det(invT))
	val <- lpsi - loginvT
	return(val)
	}

invSigma <- function(tlam=NULL,lam=NULL, psi=NULL){
	if ( is.null(dim(lam)) ) q <- 1
	else q <- ncol(lam)
	invP <- diag(1/psi, length(psi)) 
	Q <- diag(rep(1,q)) + tlam %*% invP %*% lam
	R <- invP %*% lam
	val <- invP - R %*% chol2inv(chol(Q)) %*% t(R)
	return(val)
	}
	
###################
###EM COMPONENTS###
###################

cllik <- function(ss, gpar){
	zlog <- matrix(0, nrow=nrow(ss[[1]]$y), ncol=length(gpar$pi))
	for (k in 1:length(gpar$pi) ) {
		tlam = t(gpar[[k]]$lambda)
		invS <- invSigma(tlam=tlam,lam=gpar[[k]]$lambda,psi=gpar[[k]]$psi)
		logdetS <- logdetSigma(tlam=tlam,lam=gpar[[k]]$lambda, psi=gpar[[k]]$psi)
		logdetSS <- apply(ss[[k]]$SS, 3, function(z) {
			k <- diag(z) !=0	
			if (sum(k)== 1) val = log(z[k,k])
			else val = log(det(z[k,k]))
			return(val)	} )	
		zlog[,k] <- - (logdetS-logdetSS)/2 -  mahalanobis(x=ss[[k]]$y, center=gpar[[k]]$mu, cov= invS, inverted=TRUE)/2
		}
	z <- exp(zlog)
	if ( length(gpar$pi)  >1 ) w = apply(sweep(z,2,gpar$pi, "*"),1,sum)
	else w = z	
	val = sum(log(w), na.rm=TRUE)
	return(val)
	}

getall <- function(loglik){
	if (length(loglik) >= 3) {
		n <- length(loglik)
		lm1 <- loglik[n]
		lm  <- loglik[(n-1)]
		lm_1  <- loglik[(n-2)]
		am <- (lm1 - lm)/(lm - lm_1)
		lm1.Inf <- lm + (lm1 - lm)/(1-am)
		val <- lm1.Inf - lm	
		if ( !is.na(val) ){ 
			val <- abs(val) 
			} else {
				val <- 0
				}
		} else {
			val <- 1		
			}
	return( val )
	}

getss <- function(data, gpar){
	ss <- list()
	for (i in 1:length(gpar$pi)) ss[[i]] <- initial.ySg(data, par= gpar[[i]], weights=NULL)
	return(ss)
	}

MAP <- function(data, gpar, ss=NULL){
	if (is.null(ss)) ss <- getss(data, gpar)
	w <- weights(data=data, gpar=gpar, v=1, ss=ss)
	z <- apply(w, 1, function(z) { (1:length(z))[z==max(z)] })
	return(z)	
	}

BIC <- function(like,n,g,p,q){
	val <- list()
	c <- like[!is.na(like)]
	val$c1 <- length(c)
	val$v1 <- like[val$c1]
	val$v2 <- g-1 + g*(p + p) + p*q - q*(q-1)/2
	val$bic <- 2*val$v1 - val$v2*log(n)
	return(val)
	}
	
ICL <- function(bic=NULL,weights=NULL){
	constraint <- sum(log(apply(weights,1,max)))
	val <- bic + constraint
	return(val)
	}

################
###EM UPDATES###
################

weights <- function(data=NULL, gpar=NULL, v=1, ss=NULL){
	zlog <- matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
	if ( length(gpar$pi) > 1 ) {
	for (k in 1:length(gpar$pi) ) {
		invS    = invSigma(tlam=gpar[[k]]$tlambda,lam=gpar[[k]]$lambda,psi=gpar[[k]]$psi)
		logdetS = logdetSigma(tlam=gpar[[k]]$tlambda,lam=gpar[[k]]$lambda, psi=gpar[[k]]$psi)	
		logdetSS = apply(ss[[k]]$SS, 3, function(z) {
			k = diag(z) !=0	
			if (sum(k)== 1) val = log(z[k,k])
			else  val= log(det(z[k,k]))
			return(val)	 } )	
		zlog[,k] = - (logdetS-logdetSS)/2 -  mahalanobis(x=ss[[k]]$y, center=gpar[[k]]$mu, cov= invS, inverted=TRUE)/2
		}
	w <- t(apply(zlog, 1, function(z,wt,v) {
			x=exp(v* (z + log(wt)) );
			return(x/sum(x) ) 
			}, wt=gpar$pi,v=1 ))
		} else {
			w = matrix(1, nrow=nrow(data),ncol=1)		
			}					
	return(w)
	}

initial.ySg <- function(x, par=NULL, weights=NULL, invS=NULL, ss=NULL) {
	if (is.null(weights)) weights=rep(1,nrow(x))
	if (is.null(invS)) invS = invSigma(tlam=par$tlambda, lam=par$lambda, psi=par$psi)
	q = ncol(par$lam)
	xnew <- t(apply(x,1,function(z, mu=NULL, lam=NULL, tlam=NULL, psi=NULL) {
			missing = is.na(z)		
			x1 = z[!missing]			
			mu1 = mu[!missing]
			mu2 = mu[missing]
			lam1 = matrix(lam[!missing,], nrow=sum(!missing), ncol=q)
			tlam1 = matrix(tlam[,!missing],nrow=q,ncol=sum(!missing))
			lam2 = matrix(lam[missing,], nrow=sum(missing), ncol=q)
			psi1 = psi[!missing]
			psi2 = psi[missing]
			invS11 = invSigma(tlam=tlam1,lam=lam1,psi=psi1)
			x2 = mu2 + lam2 %*% tlam1 %*% invS11 %*% (x1 - mu1)		
			x = z
			x[missing] = as.numeric(x2)
			return(x)	
			}, mu=par$mu, lam=par$lam, tlam=par$tlam, psi=par$psi)
			)
	n <- nrow(x)
	R <- matrix(0, nrow=ncol(x), ncol=ncol(x))
	RR =  array(0, dim=c(ncol(x), ncol(x), nrow(x)))
	mu <- par$mu
	tlam <- par$tlambda
	lam <- par$lambda
	psi <- par$psi
	sigma <- lam %*% tlam + diag(psi, length(psi))
	for (i in 1:n){
		missing = is.na(x[i,])
		Rt = sigma[missing, missing] - sigma[missing, !missing] %*% solve(sigma[!missing, !missing]) %*% sigma[!missing,missing]
		RR[missing,missing,i] = Rt
		R = R + RR[,,i]* weights[i]
		}	
	sumw = sum(weights)
	R = R/sumw
	if (is.null(ss)) ss = list()
	ss$ybar = apply(xnew, 2, weighted.mean, w=weights )
	ss$y = xnew
	ss$S = R + cov.wt(xnew, center=mu, wt=weights, method="ML" )$cov
	ss$SS = RR
	return(ss)
	}

update.ySg <- function(x, par=NULL, weights=NULL, ss=NULL){
	if (is.null(weights)) weights=rep(1,nrow(x))
	if (is.null(par$invS)) par$invS = invSigma(lam=par$lambda, psi=par$psi)
	y <- ss$y
	S <- ss$S
	SS <- ss$SS
	sigma <- par$lambda %*% t(par$lambda) + diag(par$psi, length(par$psi))	
	mu <- par$mu	
	xis.na <- is.na(x)
	cF <- matrix(0, ncol(x), ncol(x) )
	invS <- par$invS
	for(j in 1:ncol(x)){
		missing <- is.na(x[,j])
		if (sum(missing) > 0 ) {
		Rj <- -1*invS[j,-j]/invS[j,j]
		y[missing, j] <- mu[j] + Rj %*% t(sweep(as.matrix(y[missing,-j]), 2, mu[-j]))
		aj <- as.numeric(sum( invS[j,-j] * sigma[-j,j])  )
		if (  aj == 1  ) invSj = ginv(sigma[-j,-j])
		else  invSj = ( diag(1,ncol(x)-1) + outer(invS[j,-j], sigma[-j,j]) /as.numeric( 1 - aj )) %*% invS[-j,-j]
			for (k in 1:nrow(x) ) {
				if ( missing[k] ) {
					zerok =  xis.na[k,]			
					val2 = sigma[j, zerok] - t(sigma[-j,j]- SS[-j,j,k]) %*% invSj %*% (sigma[-j, zerok]- SS[-j, zerok,k]) 
					SS[j, zerok,k] = val2
					SS[zerok,j,k] = SS[j, zerok,k]
					}
				}
			}
		}
	cF = matrix(0, ncol(x), ncol(x) )
	for (i in 1:nrow(x)) cF = cF + SS[,,i]* weights[i]
	cF = cF/sum(weights)	
	mu <- apply(y, 2, weighted.mean, w=weights )
	Snew <- cov.wt(y, center=mu, wt=weights, method="ML" )$cov + cF
	ss <- list()
	ss$ybar <- mu
	ss$y <- y
	ss$S <- Snew
	ss$SS <- SS
	return(ss)
	}
                                                                                                                                                                                                                                                                                                                                                                                                                                                           
getmatRT <- function(x=NULL, par=NULL, weights=NULL, ssk=NULL){
	if (is.null(weights)) weights <- rep(1,nrow(x))
	mu <- par$mu
	lam <- par$lambda
	tlam <- par$tlambda
	psi <- par$psi
	n <- nrow(x)
	sumw <- sum(weights)
	R <- ssk$S
	Sigma <- lam%*%tlam + diag(psi)
	invS <- par$invS
	if (is.null(invS)) invS <- invSigma(tlam=tlam, lam=lam, psi=psi)
	q <- ncol(lam)
	beta <- tlam%*%invS
	BR <- matrix(beta%*%R, ncol=nrow(par$lambda), nrow=ncol(par$lambda))
	T <- diag(rep(1,q)) - beta%*%lam + BR%*%t(beta)
	newRT <- list(R=R, T=T, BR=BR, ng=sumw)
	return(newRT)
	}	

EMgrLstep <- function(data=NULL, gpar=NULL, ss=NULL, label=NULL) {
	w <- weights(data=data, gpar=gpar, ss=ss)
	if (!is.null(label)) w <- combinewk(weights=w, label= label)
	p <- nrow(gpar[[1]]$lambda); q <- ncol(gpar[[1]]$lambda); G <- length(gpar$pi);
	matRTB <- list()
	PT <- Diagonal(n=p*q, x=rep(0,p*q))
	PBR <- matrix(0, ncol=q, nrow=p);
	for(k in 1:G){
		ss[[k]] <- update.ySg(x=data, par=gpar[[k]], weights=w[,k], ss=ss[[k]] )
		gpar[[k]]$mu <- ss[[k]]$ybar
		matRTB[[k]] <- getmatRT(x=data, par=gpar[[k]], weights=w[,k], ssk=ss[[k]])
		wg <- matRTB[[k]]$ng/nrow(data)
		PT <- PT + wg*kronecker( Diagonal(n=p, x=1/gpar[[k]]$psi), matRTB[[k]]$T)
		PBR <- PBR + wg*( diag(1/gpar[[k]]$psi) %*% t(matRTB[[k]]$BR) )
		}
	lambda <- as.matrix(chol2inv(chol(forceSymmetric(PT))))%*%as.numeric(t(PBR))
	lambda <- t(matrix(lambda, nrow=q, ncol=p))
	tlambda <- t(lambda)
	for(k in 1:G){
		gpar[[k]]$lambda <- lambda
		gpar[[k]]$tlambda <- tlambda
		gpar[[k]]$psi <- diag( matRTB[[k]]$R - 2*(lambda %*% matRTB[[k]]$BR) + lambda%*%matRTB[[k]]$T%*%tlambda ) 
		gpar[[k]]$invS <- invSigma(tlam=tlambda,lam=lambda, psi=gpar[[k]]$psi)		
		}	
	gpar$pi <- apply(w,2,mean)
	val <- list(gpar=gpar, ss=ss)
	return(val)
	}  

EM <- function(data=NULL, max.iter=NULL, epsilon=NULL, G=NULL, q=NULL, print=print ) {
	loglik <- numeric(max.iter)
	all  <- numeric(max.iter) + 1 	
	gpar = igpar(data,g=G,q=q)
	w = matrix( runif( length(gpar$pi)*nrow(data) ), nrow=nrow(data), ncol=length(gpar$pi))
	w = t(apply(w,1,function(z) {z/sum(z)}))
	ss = list()
	if(G > 1) for(g in 1:length(gpar$pi)) ss[[g]] = initial.ySg(x=data, par=gpar[[g]], weights=w[,g]) 
	else ss[[1]] = initial.ySg(x=data, par=gpar[[1]], weights=NULL)
	i <- 1
	temp = EMgrLstep(data=data, gpar=gpar, ss=ss)
	loglik[i] = cllik(ss=temp$ss,gpar=gpar)
	gpar = temp$gpar
	ss   = temp$ss
 	if(print == TRUE) cat(paste("Iteration ",i,". The Log-Likelihood is ",loglik[i],sep=""),"\n")
	while( (all[i] > epsilon) & (i < max.iter) ){
		i <- i + 1
		temp = EMgrLstep(data=data, gpar=gpar, ss=ss)
		loglik[i] = cllik(ss=temp$ss,gpar=gpar)
 		gpar = temp$gpar
		ss   = temp$ss
		all[i]    = getall(loglik[1:i]) 
		if(print == TRUE) cat(paste("Iteration ",i,". The Log-Likelihood is ",loglik[i],sep=""),"\n")
		if(loglik[i] < loglik[i-1]) stop(paste("The",G,"component mixture with",q,"latent factors did not fit because the likelihood decreased on iteration.",i),call.=FALSE)
		}
	w <- weights(data=data,gpar=gpar,ss=ss)
	val = list(loglik=loglik[1:i], all=all[1:i], gpar=gpar, ss=ss, w=w)
	return(val)
	}

EMn <- function(data=NULL, igpar=NULL, G=NULL, q=NULL, n=10, label=label){
	loglik <- numeric(n)
	gpar <- igpar
	if (length(gpar$pi) == 1) w = NULL
	else {
		w = matrix( runif( length(gpar$pi)*nrow(data) ), nrow=nrow(data), ncol=length(gpar$pi))
		w = t(apply(w,1,function(z) {z/sum(z)}))
		}
	ss = list()
	for (g in 1:length(gpar$pi)) ss[[g]] = initial.ySg(x=data, par=gpar[[g]], weights=w[,g])
	for (i in 1:n) {			
 		temp = EMgrLstep(data=data, gpar=gpar, ss=ss, label=label) # if (i %% 100 ==0 ) 
		loglik[i] = cllik(ss=temp$ss,gpar=gpar)
		gpar = temp$gpar
		ss   = temp$ss
		}
	val = list(loglik=loglik, all=all, gpar=gpar, ss=ss)
	return(val)
	}

###############
###Algorithm###
###############

gety2 <- function(data, gpar, ss){
	w <- weights(data=data, gpar=gpar, ss=ss)
	y <- matrix(0,nrow=nrow(data), ncol=ncol(data))
	for (i in 1:length(gpar$pi)) y = y + sweep(ss[[i]]$y,1,w[,i],"*")
	return(y)
	}

getu2 <- function(data, gpar, ss){
	out <- list()
	w <- weights(data=data, gpar=gpar, ss=ss)
	u.g <- u.zg <- matrix(0, nrow=nrow(data), ncol=length(gpar$pi) ) 
	u <- matrix(0, nrow=nrow(data), ncol=ncol(gpar[[1]]$lambda) )
	for (i in 1:length(gpar$pi)){ 
		u.g[,i] <- (ss[[i]]$y-gpar[[i]]$mu) %*% gpar[[i]]$invS %*% (gpar[[i]]$lambda)
		u.zg[,i] <- sweep(ss[[i]]$y-gpar[[i]]$mu,1,w[,i],"*") %*% gpar[[i]]$invS %*% (gpar[[i]]$lambda) # Change 1
		u = u + sweep(ss[[i]]$y-gpar[[i]]$mu,1,w[,i],"*") %*% gpar[[i]]$invS %*% (gpar[[i]]$lambda)  # Change 2
		}
	out$u.g <- u.g
	out$u.zg <- u.zg
	out$u <- u
	return(out)
	}

summary.NA.fact <- function(gpar) {
	p <- nrow(gpar[[1]]$lambda)
	G <- length(gpar$pi)
	mu <- matrix(0, nrow=p,ncol=G)
	psi <- matrix(0, nrow=p,ncol=G)
	for(g in 1:G){
		mu[,g] <- gpar[[g]]$mu
		psi[,g] <- gpar[[g]]$psi
		}
	lambda <- gpar[[1]]$lambda
	val <- list(mu=mu, psi=psi, lambda=lambda, pi=gpar$pi)
	return(val)		
	}

CUUimpute <-
function(x, G=1:3, q=1:2, epsilon=1e-2, max.iter=10000, known=NULL, print=TRUE) {
	n <- nrow(x); p <- ncol(x)
	mG <- max(G); mq <- max(q)
	icl <- bic <- matrix(NA,mG,mq); 
	colnames(bic) <- paste("q=",1:mq,sep="")
	rownames(bic) <- paste("G=",1:mG,sep="")
	iclresult <- fit <- list()
	for(i in sort(G)){
		fit[[i]] <- list()
		for(j in sort(q)){
			if(print == TRUE){
				if(j == 1) cat(paste("Fitting a",i,"component mixture model with",j,"latent factor","\n"))
				else cat(paste("Fitting a",i,"component mixture model with",j,"latent factors","\n"))
				}
			fit[[i]][[j]] <- try( EM(data=x, epsilon=epsilon, max.iter=max.iter, G=i, q=j, print=print) )
			if( is.list(fit[[i]][[j]]) ){
				l.like <- length(fit[[i]][[j]]$loglik)
				if(l.like == max.iter){
					if(j == 1) cat(paste("The",i,"component mixture model with",j,"latent factor used all",max.iter,"iterations.","\n"))
					else cat(paste("The",i,"component mixture model with",j,"latent factors used all",max.iter,"iterations.","\n"))
					}
				calcs <- BIC(like=fit[[i]][[j]]$loglik,n=n,g=i,p=p,q=j)
				bic[i,j] <- calcs$bic
				icl[i,j] <- ICL(bic=bic[i,j],weights=fit[[i]][[j]]$w)
				}	
			}
		}
	if(sum(is.na(bic)) != mG*mq){
		bic <- bic[G,q]
		bic.choice <- which(max(bic,na.rm=TRUE) == bic,arr.ind=TRUE)
		if( length(bic) == 1 ) { g.bic <- G; q.bic <- q; final.bic <- bic }
		if( length(bic) > 1 && length(q) == 1 ) { q.bic <- q; g.bic <- G[bic.choice]; final.bic <- bic[bic.choice] }
		if( length(bic) > 1 && length(G) == 1 ) { g.bic <- G; q.bic <- q[bic.choice]; final.bic <- bic[bic.choice] }
		if( length(bic.choice) > 1 ) { g.bic <- G[bic.choice[1]]; q.bic <- q[bic.choice[2]]; final.bic <- bic[bic.choice[1],bic.choice[2]] }
		bestfit <- paste("The best fitting model, as selected by the BIC, is a",g.bic,"component mixture with",q.bic,"latent factors")
		yhat <- gety2(data=x,gpar=fit[[g.bic]][[q.bic]]$gpar,ss=fit[[g.bic]][[q.bic]]$ss)
		u <- getu2(data=x,gpar=fit[[g.bic]][[q.bic]]$gpar,ss=fit[[g.bic]][[q.bic]]$ss)
		gpar <- summary.NA.fact(gpar=fit[[g.bic]][[q.bic]]$gpar)
		zig <- weights(data=x,gpar=fit[[g.bic]][[q.bic]]$gpar,ss=fit[[g.bic]][[q.bic]]$ss)
		map <- MAP(data=x,gpar=fit[[g.bic]][[q.bic]]$gpar,ss=fit[[g.bic]][[q.bic]]$ss)
		if(sum(is.na(icl)) != mG*mq){
			icl <- icl[G,q]	
			icl.choice <- which(max(icl,na.rm=TRUE) == icl,arr.ind=TRUE)
			if( length(icl) == 1 ) { g.icl <- G; q.icl <- q; final.icl <- icl }
			if( length(icl) > 1 && length(q) == 1 ) { q.icl <- q; g.icl <- G[icl.choice]; final.icl <- icl[icl.choice] }
			if( length(icl) > 1 && length(G) == 1 ) { g.icl <- G; q.icl <- q[icl.choice]; final.icl <- icl[icl.choice] }
			if( length(icl.choice) > 1 ) { g.icl <- G[icl.choice[1]]; q.icl <- q[icl.choice[2]]; final.icl <- icl[icl.choice[1],icl.choice[2]] }
			iclresult$bestfit <- paste("The best fitting model, as selected by the ICL, is a",g.icl,"component mixture with",q.icl,"latent factors")
			iclresult$allicl <- icl
			iclresult$icl <- final.icl
			iclresult$yhat <- gety2(data=x,gpar=fit[[g.icl]][[q.icl]]$gpar,ss=fit[[g.icl]][[q.icl]]$ss)
			iclresult$u <- getu2(data=x,gpar=fit[[g.icl]][[q.icl]]$gpar,ss=fit[[g.icl]][[q.icl]]$ss)
			print("*******")
			iclresult$gpar <- summary.NA.fact(gpar=fit[[g.icl]][[q.icl]]$gpar)
			iclresult$zig <- weights(data=x,gpar=fit[[g.icl]][[q.icl]]$gpar,ss=fit[[g.icl]][[q.icl]]$ss)
			iclresult$map <- MAP(data=x,gpar=fit[[g.icl]][[q.icl]]$gpar,ss=fit[[g.icl]][[q.icl]]$ss)
			iclresult$G <- ncol(iclresult$zig)
			iclresult$q <- ncol(iclresult$gpar$lambda)
			}
		if( !is.null(known) ) { 
			true <- known
			class.table <- table(true,map)
			map <- iclresult$map
			iclresult$class.table <- table(true,map)
			val <- list( allbic=bic, bic=final.bic, G=ncol(zig), q=ncol(gpar$lambda), loglik=fit[[g.bic]][[q.bic]]$loglik, gpar=gpar, yhat=yhat, u=u$u, zig=zig, map=map, class.table=class.table, iclresult=iclresult )
			}
		else val <- list( allbic=bic, bic=final.bic, G=ncol(zig), q=ncol(gpar$lambda), loglik=fit[[g.bic]][[q.bic]]$loglik, gpar=gpar, yhat=yhat, u=u$u, zig=zig, map=map, iclresult=iclresult )
		print("********")
		} else stop("It appears that not one of the mixtures was able to fit this data",call.=FALSE)
	class(val) <- "CUUimpute"
	val
	}

print.CUUimpute <-function(x, ...){
	cat("Results","\n\n")
	cat("Log-Likelihood: ",x$loglik[length(x$loglik)],"\n\n")
	cat(paste("According to the BIC (",round(x$bic,digits=3),")",sep=""),"\n")
	cat(paste("# of components:            ",length(x$gpar$pi),"\n"))
	cat(paste("# of latent factors:        ",ncol(x$gpar$lambda),"\n"))
	cat(paste("According to the ICL (",round(x$iclresult$icl,digits=3),")",sep=""),"\n")
	cat(paste("# of components:            ",length(x$iclresult$gpar$pi),"\n"))
	cat(paste("# of latent factors:        ",ncol(x$iclresult$gpar$lambda),"\n"))
	}

summary.CUUimpute <- function(object, ...){
	x <- object
	cat("Results","\n\n")
	cat("Log-Likelihood: ",x$loglik[length(x$loglik)],"\n\n")
	cat(paste("According to the BIC (",round(x$bic,digits=3),")",sep=""),"\n")
	cat(paste("# of components:            ",length(x$gpar$pi),"\n"))
	cat(paste("# of latent factors:        ",ncol(x$gpar$lambda),"\n"))
	cat(paste("According to the ICL (",round(x$iclresult$icl,digits=3),")",sep=""),"\n")
	cat(paste("# of components:            ",length(x$iclresult$gpar$pi),"\n"))
	cat(paste("# of latent factors:        ",ncol(x$iclresult$gpar$lambda),"\n"))
	}