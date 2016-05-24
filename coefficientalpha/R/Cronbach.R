## The package call the rsem package for robust analysis

## R function for Cronbach alpha
alpha<-function(y, varphi=0.1, se=FALSE, test=TRUE, complete=FALSE, auxiliary=NULL, drop, silent=TRUE){
	## get the dim of y
	p<-ncol(y)
	
	## combine data
	if (!is.null(auxiliary)){
		if (!(is.matrix(auxiliary) | is.data.frame(auxiliary))) stop('auxiliary should be an matrix or a data frame')
		y.aux<-cbind(y, auxiliary)
		p.aux<-ncol(auxiliary)
	}else{
		y.aux<-y
	}
	
	## remove some cases
	if (!missing(drop)){
		drop <- -drop
		obsind <- 1:nrow(y.aux)
		obsind <- obsind[drop]
		y.aux<-y.aux[obsind, ]
	}
	if (complete) y.aux<-na.omit(y.aux)	

	n<-nrow(y.aux)
	rownames(y.aux)<-1:n
	
	pat<-rsem.pattern(y.aux)
	cov1<-rsem.emmusig(pat,varphi=varphi)
	
	weight<-rep(1, n)
	if (varphi>0) weight<-cov1$weight
	
	## covariance for y only
	sigma<-cov1$sigma[1:p, 1:p]
	mu<-cov1$mu[1:p]

	alpha<-p/(p-1)*(1-sum(diag(sigma))/sum(sigma))

	## robust standard error
	
	se.alpha<-NA
	if (se | test){
		s<-rsem.Ascov(pat, cov1, varphi=varphi)	
		## get out the Gamma matrix
		vname<-rownames(sigma)
		gname<-NULL				
		for (i in 1:p){
			for (j in i:p){
				gname<-c(gname, paste(vname[i], ".", vname[j], sep=""))
			}
		}	
		gmname<-c(vname, gname)
	}
	if (se){	
		firstd<-NULL
		
		for (i in 1:p){
			for (j in i:p){
				if (i==j){
					temp<--p/(p-1)*(sum(sigma)-sum(diag(sigma)))/((sum(sigma))^2)
				}else{
					temp<- 2*p/(p-1)*sum(diag(sigma))/((sum(sigma))^2)
				}
				firstd<-c(firstd, temp)
			}
		}
		
		## take out the part of gamma for the covariance of y
		if (is.null(auxiliary)){
			q<-p+p*(p+1)/2
			gamma<-s$Gamma[(p+1):q, (p+1):q]
		}else{
			gamma<-s$Gamma[gname, gname]
		}
		
		se.alpha<-sqrt(firstd%*%gamma%*%t(t(firstd)))/sqrt(n)
	}
	
	## conduct the test of tau-equivalent
	fit<-NULL
	if (test){
		if (is.null(auxiliary)){
			gamma<-s$Gamma
		}else{
			gamma<-s$Gamma[gmname, gmname]
		}

		vname<-rownames(sigma)
		## model for tau-equivalent
		model<-paste('f =~ 1*', vname[1])
		for (i in 2:length(vname)){
			model<-paste(model, '+ 1*', vname[i])
		}
	
		## estimate the model using rsem
		cfa.res<-cfa(model, sample.cov=sigma, sample.mean=mu, meanstructure = TRUE, sample.nobs=n)

		robust.fit <- rsem.fit(cfa.res, gamma, cov1)
		fit<-robust.fit[[4]][[1]]
		if (silent){
			cat("Test of tau-equavilence\n")
			cat("The robust F statistic is ", round(fit[1],3), "\n")
			cat("with a p-value ", round(fit[4],4), "\n")
			if (fit[4]<.05) cat("**The F test rejected tau-equavilence**\n\n")
		}
	}
	
	if (varphi>0) prop.down<-(n-sum(weight$w1==1))/n*100
	
	if (silent) cat('The alpha is ', alpha, sep='')
	if (se & silent) cat(' with the standard error ', se.alpha, sep='')
	if (silent) cat('.\n')
	
	if (varphi>0 & silent) cat('About ', round(prop.down,2), '% of cases were downweighted.\n', sep='')
	
	alphaobject<-list(alpha=alpha, se=se.alpha, weight=weight, y=y, auxiliary=auxiliary, varphi=varphi, musig=cov1, fit=fit)
	class(alphaobject)<-'alpha'
	invisible(alphaobject)	
}


plot.alpha<-function(x, type="weight", profile=5, interval=0.01, center=TRUE, scale=FALSE, w1=FALSE, numbered=FALSE, pos='topright', ...){
	## type: weight, profile, diagnosis
	res<-x
	y<-res$y
	outcase.temp<-sum(res$weight$w1<1)
	## find the smallest weight
	if (profile==0){
		outcase<-5
		if (outcase.temp<5) outcase<-outcase.temp
	}
	if (outcase.temp<5) profile<-outcase.temp
		
	## find the smallest cases
	temp<-sort(res$weight$w1)
	
	if (profile==0){
		idx<-which(res$weight$w1 <= temp[outcase])
	}else{
		idx<-which(res$weight$w1 <= temp[profile])
	}
	
	if (substr(type,1,1)=="w"){
		if (w1){
			par(mfrow=c(2,1))
			plot(res$weight$w1, ylim=c(0,1), xlab='Case number', ylab='Weight w1', main='Weights for mean estimates',...)
			text(idx, res$weight$w1[idx], idx, pos=1, col='red', ...)
			plot(res$weight$w2, ylim=c(0,max(res$weight$w2)), xlab='Case number', ylab='Weight w2', main='Weights for covariance estimates',...)		
			text(idx, res$weight$w2[idx], idx, pos=1, col='red', ...)
			par(mfrow=c(1,1))
		}else{
			plot(res$weight$w2, ylim=c(0,max(res$weight$w2)), xlab='Case number', ylab='Weight w2',...)		
			text(idx, res$weight$w2[idx], idx, pos=1, col='red', ...)
		}
	} 
	
	if (substr(type,1,1)=="p"){
		## profile plot		
		if (center){
			if (scale){
				for (i in 1:nrow(y)){
				   y[i, ]<-(y[i, ]-res$musig$mu)/sqrt(diag(res$musig$sigma))
			   }
			}else{
				for (i in 1:nrow(y)){
				   y[i, ]<-y[i, ]-res$musig$mu
				}
			}			
		}
		l.min<-min(y, na.rm=TRUE)
		l.max<-max(y, na.rm=TRUE)
	
		p<-ncol(y)
		
		if (center){
			plot(1:p, rep(0,p), type='l', ylim=c(l.min, l.max), lwd=3, ylab='Score', xlab='', axes=FALSE, col='black', ...)
		}else{
			plot(1:p, res$musig$mu, type='l', ylim=c(l.min, l.max), lwd=3, ylab='Score', xlab='', axes=FALSE, col='black', ...)
		}
		box()
		axis(1, at=1:p, labels=names(y), las=2, ...)
		axis(2, ...)
		
		if (!is.null(pos)){
			ltyno<-1	
			idx.type<-NULL	
			for (i in idx){
				lines(1:p, y[i, ], lty=ltyno, ...)
				if (numbered & center) text(1, y[i,1], i)
				ltyno<-ltyno+1
				## type of outliers
				temp.type<-'(O)'
				if (max(y[i, ], na.rm=TRUE)<0)  temp.type<-'(L-)'
				if (min(y[i, ], na.rm=TRUE)>0)  temp.type<-'(L+)'
				idx.type<- c(idx.type, temp.type)
			}
			legend(pos, legend=paste(idx, idx.type), lty=1:ltyno)
		}
	}
	
	if (substr(type,1,1)=="d"){
		varphi<-res$varphi
		phi<-seq(0, varphi, by=interval)
		alpha.diag<-NULL
		for (i in 1:length(phi)){
			alpha.diag<-c(alpha.diag, alpha(y, phi[i], auxiliary=res$auxiliary, test=FALSE, se=FALSE)$alpha)
		}
		plot(phi, alpha.diag, type='o', xlab='varphi', ylab='alpha', ylim=c(0,1), ...)
		points(phi, alpha.diag, bg='black', pch=21, cex=1, lwd=1)
	}
}

summary.alpha<-function(object, type='raw', prob=.95, ...){
    if (prob > .5) prob <- 1 - prob
	res<-object
	cat("The estimated alpha is \n")
	
	t0.txt <- sprintf("  %-20s", "alpha")
	t1.txt <- sprintf("  %10.3f", res$alpha)
	cat(t0.txt, t1.txt, "\n", sep="")
  
	if (!is.na(res$se)){
		t0.txt <- sprintf("  %-20s", "se")
		t1.txt <- sprintf("  %10.3f", res$se)
		cat(t0.txt, t1.txt, "\n", sep="")
  
		## output p-value, too
		z<-abs(res$alpha/res$se)
		pvalue<-(1-pnorm(z))
		t0.txt <- sprintf("  %-20s", "p-value (alpha>0)")
		t1.txt <- sprintf("  %10.3f", pvalue)
		cat(t0.txt, t1.txt, "\n", sep="")
  
		## output confidence interval
		if (type=='logit'){
			lambda<-log(res$alpha/(1-res$alpha))
			se.lambda<-res$se/(res$alpha*(1-res$alpha))
			ci<-c( 1/(1+exp(-lambda+abs(qnorm(prob/2))*se.lambda)),  1/(1+exp(-lambda-abs(qnorm(prob/2))*se.lambda)))
		}else{
			ci<-res$alpha+c(-1,1)*abs(qnorm(prob/2))*res$se
			if (ci[2]>1) ci[2]<-1.000
			if (ci[1]<0) ci[1]<-0.000
		}
  
		t0.txt <- sprintf("  %-20s", "Confidence interval")
		t1.txt <- sprintf("  %10.3f", ci)
		cat(t0.txt, t1.txt, "\n\n", sep="")
	}
	
	## output the fit statistic if existing
	if (!is.null(res$fit)){
		cat("\nTest of tau-equivalence \n")
		cat("  The robust F statistic is ", round(res$fit[1],3), "\n")
		cat("  with a p-value ", round(res$fit[4],4), "\n")
		if (res$fit[4]<.05){
			cat("  Note. The robust test rejected the tau-equivalence assumption.\n")
		}else{
			cat("  Note. The robust test failed to reject the tau-equivalence assumption.\n")
		}
	}
}

omega<-function(y, varphi=0.1, se=FALSE, test=TRUE, complete=FALSE, auxiliary=NULL, drop, silent=TRUE){
	## get the dim of y
	p<-ncol(y)
	
	## combine data
	if (!is.null(auxiliary)){
		if (!(is.matrix(auxiliary) | is.data.frame(auxiliary))) stop('auxiliary should be an matrix or a data frame')
		y.aux<-cbind(y, auxiliary)
		p.aux<-ncol(auxiliary)
	}else{
		y.aux<-y
	}
	
	## remove some cases
	if (!missing(drop)){
		drop <- -drop
		obsind <- 1:nrow(y.aux)
		obsind <- obsind[drop]
		y.aux<-y.aux[obsind, ]
	}
	if (complete) y.aux<-na.omit(y.aux)	

	n<-nrow(y.aux)
	rownames(y.aux)<-1:n
	
	pat<-rsem.pattern(y.aux)
	cov1<-rsem.emmusig(pat,varphi=varphi)
	
	weight<-rep(1, n)
	if (varphi>0) weight<-cov1$weight
	
	## covariance for y only
	sigma<-cov1$sigma[1:p, 1:p]
	mu<-cov1$mu[1:p]
	
	## fit a factor model with the covariance matrix and get parameter estimates and their covariance matrix
	vname<-rownames(sigma)
	
	## the factor model
	model<-paste('f =~ ', vname[1])
	for (i in 2:length(vname)){
		model<-paste(model, '+', vname[i])
	}
	
	## estimate the model using rsem
	cfa.res<-cfa(model, sample.cov=sigma, sample.mean=cov1$mu, meanstructure = TRUE, sample.nobs=n, std.lv=TRUE)

	## calculate omega
	cfa.load<-cfa.res@Fit@est[cfa.res@ParTable$op=="=~" & cfa.res@ParTable$free>0]
	cfa.psi<-cfa.res@Fit@est[cfa.res@ParTable$op=="~~" & cfa.res@ParTable$free>0]
	
	p1<-(sum(cfa.load))^2
	p2<- p1 + sum(cfa.psi)
	omega<- p1/p2

	## robust standard error
	se.omega<-NA
	if (se | test){
		s<-rsem.Ascov(pat, cov1, varphi=varphi)	
		## get out the Gamma matrix
		vname<-rownames(sigma)
		gname<-NULL				
		for (i in 1:p){
			for (j in i:p){
				gname<-c(gname, paste(vname[i], ".", vname[j], sep=""))
			}
		}	
		gmname<-c(vname, gname)
		Gamma<-s$Gamma[gmname, gmname]
	}

	if (se){		
		cfa.se<-rsem.se(cfa.res, Gamma)
		firstd<-rep(sum(cfa.load)*sum(cfa.psi)*2/(p2^2) , p)
		firstd<-c(firstd, rep( -p1/(p2^2),p))
		q<-p*2
		gamma<-cfa.se$vcov[[1]][1:q, 1:q]
		se.omega<-sqrt(firstd%*%gamma%*%t(t(firstd)))/sqrt(n)
	}
	fit<-NULL
	if (test){
		robust.fit <- rsem.fit(cfa.res, Gamma, cov1)
		fit<-robust.fit[[4]][[1]]
		if (silent){
			cat("Test of homogeneity\n")
			cat("The robust F statistic is ", round(fit[1],3), "\n")
			cat("with a p-value ", round(fit[4],4), "\n\n")
			if (fit[4]<.05) cat("**The F test rejected homogeneity**\n\n")
		}
	}
	
	if (varphi>0) prop.down<-(n-sum(weight$w1==1))/n*100
	
	if (silent)cat('The omega is ', omega, sep='')
	if (se & silent) cat(' with the standard error ', se.omega, sep='')
	if (silent)cat('.\n')
	
	if (varphi>0 & silent) cat('About ', round(prop.down,2), '% of cases were downweighted.\n', sep='')
	
	omegaobject<-list(omega=omega, se=se.omega, weight=weight, y=y, auxiliary=auxiliary, varphi=varphi, musig=cov1, fit=fit)
	class(omegaobject)<-'omega'
	invisible(omegaobject)	
}

summary.omega<-function(object, type='raw', prob=.95,...){
    if (prob > .5) prob <- 1 - prob
	res<-object
	cat("The estimated omega is \n")
	
	t0.txt <- sprintf("  %-20s", "omega")
	t1.txt <- sprintf("  %10.3f", res$omega)
	cat(t0.txt, t1.txt, "\n", sep="")
  
	if (!is.na(res$se)){
		t0.txt <- sprintf("  %-20s", "se")
		t1.txt <- sprintf("  %10.3f", res$se)
		cat(t0.txt, t1.txt, "\n", sep="")
  
		## output p-value, too
		z<-abs(res$omega/res$se)
		pvalue<-(1-pnorm(z))
		t0.txt <- sprintf("  %-20s", "p-value (omega>0)")
		t1.txt <- sprintf("  %10.3f", pvalue)
		cat(t0.txt, t1.txt, "\n", sep="")
  
		## output confidence interval
		if (type=='logit'){
			lambda<-log(res$omega/(1-res$omega))
			se.lambda<-res$se/(res$omega*(1-res$omega))
			ci<-c( 1/(1+exp(-lambda+abs(qnorm(prob/2))*se.lambda)),  1/(1+exp(-lambda-abs(qnorm(prob/2))*se.lambda)))
		}else{
			ci<-res$omega+c(-1,1)*abs(qnorm(prob/2))*res$se
			if (ci[2]>1) ci[2]<-1.000
			if (ci[1]<0) ci[1]<-0.000
		}
  
		t0.txt <- sprintf("  %-20s", "Confidence interval")
		t1.txt <- sprintf("  %10.3f", ci)
		cat(t0.txt, t1.txt, "\n\n", sep="")
	}
	
	if (!is.null(res$fit)){
		cat("\nTest of homogeneity\n")
		cat("  The robust F statistic is ", round(res$fit[1],3), "\n")
		cat("  with a p-value ", round(res$fit[4],4), "\n")
		if (res$fit[4]<.05){
			cat("  Note. The robust test rejected the assumption of  homogeneity.\n")
		}else{
			cat("  Note. The robust test failed to reject the assumption of  homogeneity.\n")
		}
	}
}

plot.omega<-function(x, type="weight", profile=5, interval=0.01, center=TRUE, scale=FALSE, w1=FALSE, numbered=FALSE, pos='topright', ...){
	## type: weight, profile, diagnosis
	res<-x
	y<-res$y
	outcase.temp<-sum(res$weight$w1<1)
	## find the smallest weight
	if (profile==5){
		outcase<-5
		if (outcase.temp<5) outcase<-outcase.temp
	}
	if (outcase.temp<5) profile<-outcase.temp
		
	## find the smallest cases
	temp<-sort(res$weight$w1)
	
	if (profile==5){
		idx<-which(res$weight$w1 <= temp[outcase])
	}else{
		idx<-which(res$weight$w1 <= temp[profile])
	}
	
	if (substr(type,1,1)=="w"){
		if (w1){
			par(mfrow=c(2,1))
			plot(res$weight$w1, ylim=c(0,1), xlab='Case number', ylab='Weight w1', main='Weights for mean estimates',...)
			text(idx, res$weight$w1[idx], idx, pos=1, col='red', ...)
			plot(res$weight$w2, ylim=c(0,max(res$weight$w2)), xlab='Case number', ylab='Weight w2', main='Weights for covariance estimates',...)		
			text(idx, res$weight$w2[idx], idx, pos=1, col='red', ...)
			par(mfrow=c(1,1))
		}else{
			plot(res$weight$w2, ylim=c(0,max(res$weight$w2)), xlab='Case number', ylab='Weight w2',...)		
			text(idx, res$weight$w2[idx], idx, pos=1, col='red', ...)
		}
	} 
	
	if (substr(type,1,1)=="p"){
		## profile plot		
		if (center){
			if (scale){
				for (i in 1:nrow(y)){
				   y[i, ]<-(y[i, ]-res$musig$mu)/sqrt(diag(res$musig$sigma))
			   }
			}else{
				for (i in 1:nrow(y)){
				   y[i, ]<-y[i, ]-res$musig$mu
				}
			}			
		}
		l.min<-min(y, na.rm=TRUE)
		l.max<-max(y, na.rm=TRUE)
	
		p<-ncol(y)
		
		if (center){
			plot(1:p, rep(0,p), type='l', ylim=c(l.min, l.max), lwd=3, ylab='Score', xlab='', axes=FALSE, col='black', ...)
		}else{
			plot(1:p, res$musig$mu, type='l', ylim=c(l.min, l.max), lwd=3, ylab='Score', xlab='', axes=FALSE, col='black', ...)
		}
		box()
		axis(1, at=1:p, labels=names(y), las=2, ...)
		axis(2, ...)
		
		if (!is.null(pos)){
			ltyno<-1	
			idx.type<-NULL	
			for (i in idx){
				lines(1:p, y[i, ], lty=ltyno, ...)
				if (numbered & center) text(1, y[i,1], i)
				ltyno<-ltyno+1
				## type of outliers
				temp.type<-'(O)'
				if (max(y[i, ], na.rm=TRUE)<0)  temp.type<-'(L-)'
				if (min(y[i, ], na.rm=TRUE)>0)  temp.type<-'(L+)'
				idx.type<- c(idx.type, temp.type)
			}
			legend(pos, legend=paste(idx, idx.type), lty=1:ltyno)
		}
	}
	
	if (substr(type,1,1)=="d"){
		varphi<-res$varphi
		phi<-seq(0, varphi, by=interval)
		omega.diag<-NULL
		for (i in 1:length(phi)){
			omega.diag<-c(omega.diag, omega(y, phi[i], auxiliary=res$auxiliary, test=FALSE, se=FALSE)$omega)
		}
		plot(phi, omega.diag, xlab='varphi', type='o', lwd=1, ylab='omega', ylim=c(0,1), ...)
		points(phi, omega.diag, bg='black', pch=21, cex=1, lwd=1)
	}
}


tau.test<-function(y, varphi=0.1, complete=FALSE, drop){
	if (!missing(drop)){
		drop <- -drop
		obsind <- 1:nrow(y)
		obsind <- obsind[drop]
		y<-y[obsind, ]
	}
	if (complete) y<-na.omit(y)	
	p<-ncol(y)
	n<-nrow(y)
	rownames(y)<-1:n
	
	pat<-rsem.pattern(y)
	cov1<-rsem.emmusig(pat,varphi=varphi)	
	ascov<-rsem.Ascov(pat, cov1, varphi=varphi)
	
	sigma<-cov1$sigma
	
	## fit a factor model with the covariance matrix and get parameter estimates and their covariance matrix
	vname<-rownames(sigma)
	
	## estimate two models
	## model for tau-equivalent
	model<-paste('f =~ 1*', vname[1])
	for (i in 2:length(vname)){
		model<-paste(model, '+ 1*', vname[i])
	}
	
	## estimate the model using rsem
	cfa.res1<-cfa(model, sample.cov=sigma, sample.mean=cov1$mu, meanstructure = TRUE, sample.nobs=n)
	
	robust.fit1 <- rsem.fit(cfa.res1, ascov$Gamma, cov1)
	fit1<-robust.fit1[[4]][[1]]
	## model with the single factor
	model<-paste('f =~ ', vname[1])
	for (i in 2:length(vname)){
		model<-paste(model, '+', vname[i])
	}
	
	## estimate the model using rsem
	cfa.res2<-cfa(model, sample.cov=sigma, sample.mean=cov1$mu, meanstructure = TRUE, sample.nobs=n, std.lv=TRUE)
	
	robust.fit2 <- rsem.fit(cfa.res2, ascov$Gamma, cov1)
	fit2<-robust.fit2[[4]][[1]]
	
	cat("Test of tau equivalent\n")
	cat("The robust F statistic is ", round(fit1[1],3), "\n")
	cat("with a p-value ", round(fit1[4],4), "\n\n")

	
	cat("Test of homogeneous items\n")
	cat("The robust F statistic is ", round(fit2[1],3), "\n")
	cat("with a p-value ", round(fit2[4],4), "\n")

	test<-list(tau=cfa.res1,homo=cfa.res2,tau.fit=fit1,homo.fit=fit2)
	class(test)<-'test'
	invisible(test)	
}

bootstrap<-function(y, type="omega", alpha=.95, nboot=1000, ci="bc", plot=FALSE, varphi=0, complete=FALSE, auxiliary=NULL, silent=FALSE){
	
	## save the results to a vector
	boot.res<-rep(NA, nboot)
	n<-nrow(y)
	if (type=="omega"){
		## for omega
		temp<-try(omega(y,varphi=varphi, se=FALSE, test=FALSE, silent=silent))
		orig.est<-temp$omega
		for (i in 1:nboot){
			bdata<-y[sample(1:n, n, replace=TRUE), ]
			temp<-try(omega(bdata,varphi=varphi, se=FALSE, test=FALSE, silent=silent))
			if (class(temp) != "try-error") {
				boot.res[i] <- temp$omega
			}
		}
	}else{
		## for alpha
		temp<-try(alpha(y,varphi=varphi, se=FALSE, test=FALSE, silent=silent))
		orig.est<-temp$alpha
		for (i in 1:nboot){
			bdata<-y[sample(1:n, n, replace=TRUE), ]
			temp<-try(alpha(bdata,varphi=varphi, se=FALSE, test=FALSE, silent=silent))
			if (class(temp) != "try-error") {
				boot.res[i] <- temp$alpha
			}
		}
	}
	
	boot.se<-sd(boot.res, na.rm=TRUE)
	if (alpha>.5) alpha<-1-alpha
	prob<-c(alpha/2, 1-alpha/2)
	if (ci=="bc"){
		z0<-qnorm(mean(boot.res<orig.est, na.rm=TRUE))
		prob<-c(pnorm(2*z0+qnorm(alpha/2)), pnorm(2*z0+qnorm(1-alpha/2)))
		boot.ci<-quantile(boot.res, prob=prob)
	}else{
		boot.ci<-quantile(boot.res, prob=prob)
	}
	
	## output the results
	if (type=="omega"){
		cat('The estimated omega is ', orig.est, "\n")
	}else{
		cat('The estimated alpha is ', orig.est, "\n")
	}
	cat('Its bootstrap se is ', round(boot.se, 3), "\n")
	cat('Its bootstrap confidence interval is [', round(boot.ci[1], 3), " , ", round(boot.ci[2], 3),"]\n")
	
	if (plot) plot(density(boot.res), xlab='Bootstrap distribution', ylab='density', main="")
	
	bootobj<-list(estimate=orig.est, se=boot.se, ci=boot.ci, boot.est=boot.res, type=type)
}
