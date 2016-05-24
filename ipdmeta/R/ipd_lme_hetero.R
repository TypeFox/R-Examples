ipd_lme_fixef.hetero <- function(n,y,sigma2,D){
	# LIST OF AGGREGATE DATA	
	qrs <- function(N,Y,sigma2){
		quadratic.forms(N,Y,sigma2,D)	
	}
	
	components <- mapply(qrs,N=n,Y=y,sigma2=sigma2,SIMPLIFY=FALSE)
	
	XWYs <- sapply(components,function(x) x$XWY)
	XWXs <- sapply(components,function(x) x$XWX)
	
	XWY <- rowSums(XWYs)
	XWX <- rowSums(XWXs)
	p <- sqrt(length(XWX))
	XWX <- matrix(XWX,p,p)
	
	V <- solve(XWX)
	beta <- V%*%XWY

	labels <- ipd_lme_labels(n)
	row.names(beta) <- labels
	row.names(V) <- labels
	colnames(V) <- labels
	
list(
	fixef = beta,
	vcov.fixef = V
)

}

ipd_lme_ranef.hetero <- function(n,y,fixef,sigma2,D){
	
	mapply(ipd_lme_ranef_study,n=n,y.bar=y,sigma2=sigma2,
				  MoreArgs=list(fixef=fixef,D=D),SIMPLIFY=FALSE)
}

ipd_lme_labels <- function(N){
	
	labels <- c("Intercept","trt")
	x.labels <- row.names(N[[1]])[-nrow(N[[1]])]
	interaction.labels <- paste("trt",x.labels,sep=":")

c(labels,x.labels,interaction.labels)
}


ipd_lme_sigma2.hetero <- function(n,y,beta,alphas,s2){

	if(is.list(s2)){
		Qs <- mapply(ipd_lme_subgroup_sigma2_study,n=n,y=y,alpha=alphas,s2=s2,MoreArgs=list(beta=beta))
	}
	else{
		Qs <- mapply(ipd_lme_subgroup_sigma2_study,n=n,y=y,alpha=alphas,s2=s2,MoreArgs=list(beta=beta))
	}
	
	N <- sapply(n,function(x)sum(c(x$trt,x$ctrl)))
	
Qs/N
}


initialize_vcov.hetero <- function(n, y, s2){
	
	# FOR RESIDUAL TAKE MEAN OF S2
	NT <- sapply(n,function(x)sum(x$trt))
	NC <- sapply(n,function(x)sum(x$ctrl))
	
	YT <- mapply(function(n,y){sum(n$trt*y$trt)/sum(n$trt)},n=n,y=y)
	YC <- mapply(function(n,y){sum(n$ctrl*y$ctrl)/sum(n$ctrl)},n=n,y=y)
	
	if(is.list(s2)){
		s2 <- mapply(function(n,y){sum((unlist(n)-1)*unlist(y))/(sum(unlist(n))-1)},n=n,y=s2)
	}
	s2[is.na(s2)] <- 0
	sigma2 <- s2
	y.bar.ctrl <-  sum((NC-1)*(YC-mean(YC))^2)/sum(NC)
	y.bar.trt <- sum((NT-1)*(((YT-mean(YT))-(YC-mean(YC))))^2)/sum(NT)
	
	list(
		sigma2 = sigma2,
		D = diag(c(y.bar.ctrl,y.bar.trt))
		)	
}

ipdlme.hetero <- function(n, y, s2,max.iter=100,tol=1e-10){
	
	vcov.init <- initialize_vcov.hetero(n, y, s2)
	beta.trace <- NULL
	alpha.trace <- NULL
	sigma.trace <- NULL
	D.trace <- NULL
	convergence.beta <- tol
	convergence.alpha <- tol
			
	# UPDATE FIXED POP AND STUDY EFFECTS
	for(j in 1:max.iter){

		beta <- ipd_lme_fixef.hetero(n,y,sigma2=vcov.init$sigma2,D=vcov.init$D)
		alpha <- ipd_lme_ranef.hetero(n,y,beta,sigma2=vcov.init$sigma2,D=vcov.init$D)
	
    	# UPDATE VARIANCE COMPONENTS
		alphas <- lapply(alpha,function(x)x$ranef)
		vcov.init$sigma2 <- ipd_lme_sigma2.hetero(n,y,s2,alphas=alphas,beta=beta$fixef)
		vcov.init$D <- ipd_lme_D(alphas)

		beta.trace <- cbind(beta.trace,beta$fixef)
		alpha.trace <- cbind(alpha.trace,unlist(alphas))
		sigma.trace <- cbind(sigma.trace,vcov.init$sigma2)
		D.trace <- cbind(D.trace,as.numeric(vcov.init$D))	
		
		
		if(j>1){
			beta.change <- ((beta.trace[,ncol(sigma.trace)-1]-beta.trace[,ncol(sigma.trace)])^2)/beta.trace[,ncol(sigma.trace)-1]^2
			alpha.change <- (alpha.trace[,ncol(sigma.trace)-1]-alpha.trace[,ncol(sigma.trace)])^2/alpha.trace[,ncol(sigma.trace)-1]^2
		    convergence.beta <- max(beta.change)
		    convergence.alpha <- max(alpha.change)
		}
		
			test <- convergence.beta<tol&convergence.alpha<tol
			if(test) break
		
	}
	
	alphas.int <- sapply(alpha,function(x)x$ranef[1])
	alphas.trt <- sapply(alpha,function(x)x$ranef[2])
	
	alphas <- cbind(alphas.int,alphas.trt)
	colnames(alphas) <- c("Intercept","Treatment")
	row.names(alphas) <- names(alpha)
	vcov.alphas <- lapply(alpha,function(x)x$vcov)
	
	colnames(vcov.init$D) <- c("Intercept","Treatment")
	row.names(vcov.init$D) <- c("Intercept","Treatment")
	
	df <- sum(sapply(n,function(x) sum(unlist(x))))-2*(nrow(alphas)-1)
	
	new("ipdlme",
			fixef = beta$fixef,
			ranef = alphas,
			vcov.fixef = beta$vcov,
			vcov.ranef = vcov.alphas,
			sigma2 = vcov.init$sigma2,
			VarCorr = vcov.init$D,
			convergence.trace =list(
				fixef = beta.trace,
				ranef = alpha.trace,
				sigma2 = sigma.trace,
				VarCorr = D.trace
			),
			converged = j<max.iter,
			n.iter = j,
			max.iter = max.iter,
			tol = tol,
			df = df
		)
}

