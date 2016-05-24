LiuParams<-function(a, b, c) {
    re <- list(df = b^2/c, muQ = a, sigmaQ = sqrt(2 * b))
    return(re)
}

pKuonen<-function (x, lambda, delta = rep(0, length(lambda))) 
{ # modified from Thomas Lumley's saddle function in survey package to include noncentrality parameters, which is licensed under GPL-2 | GPL-3
    delta <- delta[lambda > 0]
    lambda <- lambda[lambda > 0]
    if (x <= 0) 
        return(1)
    if (length(lambda) == 1) {
        pchisq(x/lambda, df = 1, ncp = delta, lower.tail = FALSE)
    }
    d <- max(lambda)
    lambda <- lambda/d
    x <- x/d
    k0 <- function(zeta) {
        -sum(log(1 - 2 * zeta * lambda))/2 + sum((delta * lambda * 
            zeta)/(1 - 2 * zeta * lambda))
    }
    kprime0 <- function(zeta) {
        sapply(zeta, function(zz) {
            sum(lambda/(1 - 2 * zz * lambda)) + sum((delta * 
                lambda)/(1 - 2 * zz * lambda) + 2 * (delta * 
                zz * lambda^2)/(1 - 2 * zz * lambda)^2)
        })
    }
    kpprime0 <- function(zeta) {
        2 * sum(lambda^2/(1 - 2 * zeta * lambda)^2) + sum((4 * 
            delta * lambda^2)/(1 - 2 * zeta * lambda)^2 + 8 * 
            delta * zeta * lambda^3/(1 - 2 * zeta * lambda)^3)
    }
    n <- length(lambda)
    if (any(lambda < 0)) {
        lmin <- max(1/(2 * lambda[lambda < 0])) * 0.99999
    }
    else if (x > sum(lambda)+sum(delta*lambda)) {
        lmin <- -0.01
    }
    else {
        lmin <- -length(lambda)*max(1+delta)/(2 * x)
    }
    lmax <- min(1/(2 * lambda[lambda > 0])) * 0.99999
    hatzeta <- uniroot(function(zeta) kprime0(zeta) - x, lower = lmin, 
        upper = lmax, tol = 1e-08)$root
    w <- sign(hatzeta) * sqrt(2 * (hatzeta * x - k0(hatzeta)))
    v <- hatzeta * sqrt(kpprime0(hatzeta))
    if (abs(hatzeta) < 1e-04) 
        NA
    else pnorm(w + log(v/w)/w, lower.tail = FALSE)
}

rareGE_joint_mc<-function(r,V,G1,G2,rho,B=10000) {
	Q<-rho*sum((t(G1)%*%r)^2)+(1-rho)*sum((t(G2)%*%r)^2)
   	lambdas <- NULL
	liuparams <- NULL
	ps <- Ts <- numeric(length(rho))
	V11<-t(G1)%*%V%*%G1
	V12<-t(G1)%*%V%*%G2
	V22<-t(G2)%*%V%*%G2
	V11V11 <- V11 %*% V11
	V11V12 <- V11 %*% V12
	V12V22 <- V12 %*% V22
	V12V12 <- V12 %*% t(V12)
	V22V22 <- V22 %*% V22
	a1 <- sum(diag(V11))
	a2 <- sum(diag(V22))
	b1 <- sum(diag(V11V11))
	b2 <- sum(diag(V12V12))
	b3 <- sum(diag(V22V22))
	c1 <- sum(diag(V11V11 %*% V11V11))
	c2 <- sum(diag(V11V12 %*% t(V11V12)))
	c3 <- sum(diag(V11V12 %*% t(V12V22)))
	c4 <- sum(diag(V12V12 %*% V12V12))
	c5 <- sum(diag(V12V22 %*% t(V12V22)))
	c6 <- sum(diag(V22V22 %*% V22V22))
	for (i in 1:length(rho)) {
		lam <- rho[i] * a1 + (1-rho[i]) * a2
		lam2 <- rho[i]^2 * b1 + 2 * rho[i] * (1-rho[i]) * b2 + (1-rho[i])^2 * b3
		lam4 <- rho[i]^4 * c1 + 4 * rho[i]^3 * (1-rho[i]) * c2 + rho[i]^2 * (1-rho[i])^2 * (4 * c3 + 2 * c4) + 4 * rho[i] * (1-rho[i])^3 * c5 + (1-rho[i])^4 * c6
		tmp.param <- LiuParams(lam, lam2, lam4)
		ps[i] <- pchisq((Q[i]-tmp.param$muQ)/tmp.param$sigmaQ*sqrt(2*tmp.param$df)+tmp.param$df, tmp.param$df, lower.tail=F)
		liuparams <- c(liuparams, list(tmp.param))
	}
	minp <- min(ps)

	for (i in 1:length(rho)) {
		Ts[i] <- (qchisq(minp, liuparams[[i]]$df, lower.tail=F)-liuparams[[i]]$df)/sqrt(2*liuparams[[i]]$df)*liuparams[[i]]$sigmaQ+liuparams[[i]]$muQ
	}

	V.cond<-V22-t(V12)%*%ginv(V11)%*%V12
	eig <- eigen(V.cond, symmetric = TRUE)
	eigval<-eig$values
	D <- diag(eigval)
	diag(D)[zapsmall(diag(D)) > 0] <- 1/sqrt(diag(D)[zapsmall(diag(D)) > 0])
	diag(D)[diag(D) <= 0] <- 0
	meanvec <- D %*% t(eig$vectors) %*% t(V12) %*% ginv(V11)

	IDX<-which(rho>=0.999)
	if(length(IDX)>0) {
		rho[IDX]<-0.999
	}

	Fcond<-function(x) {
		qtmp<-min((Ts-rho*sum(x^2))/(1-rho))
		ptmp<-pKuonen(qtmp, lambda=eigval, delta=c(meanvec%*%x)^2)
		return(ptmp)
	}

	tmpa<-mvrnorm(B, rep(0, nrow(V11)), V11)
	tmpres<-apply(tmpa, 1, Fcond)
	actualp<-mean(tmpres, na.rm=T)

	if(length(IDX)>0) {
		rho[IDX]<-1
	}

	return(list(actualp=actualp, minp=minp, rho=rho[which.min(ps)], ps=ps, info=summary(tmpres)))
}

wuweights <- function(maf) ifelse(maf>0, dbeta(maf,1,25), 0)

rareGE<-function(phenotype, genotypes, covariates, mainweights=wuweights, interweights=wuweights, family="gaussian", binomialimpute=FALSE, rho=seq(0,1,by=0.1), B=10000, INT_FIX=TRUE, INT_RAN=TRUE, JOINT=TRUE) {
	fam<-try(match.arg(family, c("gaussian", "binomial")))
	if(class(fam)=="try-error") {
		print(family)
		stop("The argument 'family' has to be either gaussian or binomial")
	}
	# phenotype is a vector, genotypes is a matrix, covariates is a matrix
        if(class(phenotype)!="numeric" && class(phenotype)!="integer") stop("phenotype should be a numeric vector!")
        n<-length(phenotype)
        if(is.data.frame(genotypes)) genotypes<-as.matrix(genotypes)
        if(!is.matrix(genotypes)) stop("genotypes should be a matrix!")
        if(nrow(genotypes)!=n) stop("Number of individuals inconsistent between phenotype and genotypes. Check your data...")
        covariates<-as.matrix(covariates)
        if(nrow(covariates)!=n) stop("Number of individuals inconsistent between phenotype and covariates. Check your data...")
        missidx<-is.na(phenotype) | apply(is.na(covariates), 1, any)
        phenotype<-phenotype[!missidx]
        genotypes<-as.matrix(genotypes[!missidx,])
        covariates<-covariates[!missidx,]
        n<-length(phenotype)
	Z<-genotypes
	nZ<-ncol(Z)
	MAF<-colMeans(Z, na.rm = TRUE)/2
	if (is.function(mainweights)) {
		weights1<-mainweights(MAF)
	} else {
                if(!is.vector(mainweights)) stop("Check the class of your variable mainweights: should be function or vector!")
                if(length(mainweights)!=nZ) stop("Number of variants inconsistent between genotypes and mainweights. Check your data...")
		weights1<-mainweights
	}
	if (is.function(interweights)) {
		weights2<-interweights(MAF)
	} else {
                if(!is.vector(interweights)) stop("Check the class of your variable interweights: should be function or vector!")
                if(length(interweights)!=nZ) stop("Number of variants inconsistent between genotypes and interweights. Check your data...")
		weights2<-interweights
	}
	# impute missing genotypes
	for(pp in 1:nZ) {
		IDX<-which(is.na(Z[,pp]))
		if(length(IDX)>0) {
			if(binomialimpute) { # random imputation based on binomial distribution
				Z[IDX,pp]<-rbinom(length(IDX), 2, MAF[pp])
			} else { # default: set missing to 0
				Z[IDX,pp]<-0
			}
		}
	}
	# now Z should not have any missing values
	G<-t(t(Z)*weights1)
	EG<-t(t(Z)*weights2)*(covariates[,1]-mean(covariates[,1]))

	# null model
	X<-cbind(rep(1,n),as.matrix(covariates))
	tmpdata<-data.frame(phenotype,id=rep(1,n),covariates,G)
	nc<-ncol(covariates)
	if(nc==1) {
		fixedmod<-"phenotype ~ covariates"
	} else {
		fixedmod<-paste("phenotype ~", paste(names(tmpdata)[3:(2+nc)],collapse=" + "))
	}
	if(nZ==1) {
		if(INT_RAN) randommod<-"~ 0 + G"
		if(INT_FIX) lmmod<-paste(fixedmod,"G",sep=" + ")
	} else {
		if(INT_RAN) randommod<-paste("~ 0 +", paste(names(tmpdata)[(nc+3):(2+nc+nZ)],collapse=" + "))
		if(INT_FIX) lmmod<-paste(fixedmod, paste(names(tmpdata)[(nc+3):(2+nc+nZ)],collapse=" + "), sep=" + ")
	}
	pINT_FIX<-pINT_RAN<-pJOINT<-pJOINTmin<-pJOINTrho<-pJOINTps<-pJOINTinfo<-NULL
	if(INT_RAN) {
		nullmod1<-glmmPQL(fixed=as.formula(fixedmod),data=tmpdata,random=list(id=pdIdent(as.formula(randommod))),family=fam,niter=300,verbose=F)
		phat<-as.numeric(predict(nullmod1, data=tmpdata, type="response", level=1))
		DELTA<-if(fam=="binomial") diag(1/phat/(1-phat)) else diag(n)
		res1<-tmpdata$phenotype - phat
		varc<-VarCorr(nullmod1)
		deltai<-if(fam=="binomial") phat*(1-phat)/as.numeric(varc[nrow(varc),1]) else rep(1/as.numeric(varc[nrow(varc),1]),n)
		Gdeltai<-G*deltai
		Di<-solve(t(G)%*%(Gdeltai)+diag(rep(1/as.numeric(varc[1,1]),nZ)))
		XSX<-t(X)%*%(X*deltai)-t(X)%*%Gdeltai%*%Di%*%(t(Gdeltai)%*%X)
		XSEG<-t(X)%*%(EG*deltai)-t(X)%*%Gdeltai%*%Di%*%(t(Gdeltai)%*%EG)
		EGSEG<-t(EG)%*%(EG*deltai)-t(EG)%*%Gdeltai%*%Di%*%(t(Gdeltai)%*%EG)
		Q1<-sum((t(EG) %*% res1 / as.numeric(varc[nrow(varc),1]))^2)
		eig1<-eigen(EGSEG-t(XSEG)%*%solve(XSX,XSEG), symmetric=T, only.values=T)
		evals1<-eig1$values[zapsmall(eig1$values)>0]
		pINT_RAN<-pchisqsum(Q1,rep(1,length(evals1)),evals1,lower.tail=F,method="sad")
	}
	if(INT_FIX) {
		X2<-cbind(X,as.matrix(G))
		nullmod2<-glm(as.formula(lmmod),data=tmpdata,family=fam)
		res2<-residuals(nullmod2, type="response")
		sigma2<-ifelse(fam=="binomial", 1, var(residuals(nullmod2, "pearson"))*(n-1)/n)
		Q2<-sum((t(EG) %*% res2)^2)/sigma2
		binvar2<-nullmod2$family$var(as.numeric(nullmod2$fitted))
		EGbEG<-t(EG*binvar2) %*% EG
		X2bEG<-t(X2*binvar2) %*% EG
		X2bX2<-t(X2*binvar2) %*% X2
		eig2<-eigen(EGbEG - t(X2bEG) %*% ginv(X2bX2) %*% X2bEG, symmetric=T, only.values=T)
		evals2<-eig2$values[zapsmall(eig2$values)>0]
		pINT_FIX<-pchisqsum(Q2,rep(1,length(evals2)),evals2,lower.tail=F,method="sad")
	}
	if(JOINT) {
		nullmod3<-glm(as.formula(fixedmod),data=tmpdata,family=fam)
		sigma2_3<-ifelse(fam=="binomial", 1, var(residuals(nullmod3, "pearson"))*(n-1)/n)
		res3<-residuals(nullmod3, type="response")/sqrt(sigma2_3)
		binvar3<-nullmod3$family$var(as.numeric(nullmod3$fitted))
		VV<-diag(binvar3)-(X*binvar3) %*% ginv(t(X*binvar3) %*% X) %*% t(X*binvar3)
		JOINT<-rareGE_joint_mc(res3, VV, G, EG, rho=rho, B=B)
		pJOINT<-JOINT$actualp
		pJOINTmin<-JOINT$minp
		pJOINTrho<-JOINT$rho
		pJOINTps<-JOINT$ps
		pJOINTinfo<-JOINT$info
	}
	return(list(pINT_FIX=pINT_FIX,pINT_RAN=pINT_RAN,pJOINT=pJOINT,pJOINTmin=pJOINTmin,pJOINTrho=pJOINTrho,pJOINTps=pJOINTps,pJOINTinfo=pJOINTinfo))
}

INT_FIX<-function(phenotype, genotypes, covariates, mainweights=wuweights, interweights=wuweights, family="gaussian", binomialimpute=FALSE) {
	fam<-try(match.arg(family, c("gaussian", "binomial")))
	if(class(fam)=="try-error") {
		print(family)
		stop("The argument 'family' has to be either gaussian or binomial")
	}
	# phenotype is a vector, genotypes is a matrix, covariates is a matrix
        if(class(phenotype)!="numeric" && class(phenotype)!="integer") stop("phenotype should be a numeric vector!")
        n<-length(phenotype)
        if(is.data.frame(genotypes)) genotypes<-as.matrix(genotypes)
        if(!is.matrix(genotypes)) stop("genotypes should be a matrix!")
        if(nrow(genotypes)!=n) stop("Number of individuals inconsistent between phenotype and genotypes. Check your data...")
        covariates<-as.matrix(covariates)
        if(nrow(covariates)!=n) stop("Number of individuals inconsistent between phenotype and covariates. Check your data...")
        missidx<-is.na(phenotype) | apply(is.na(covariates), 1, any)
        phenotype<-phenotype[!missidx]
        genotypes<-as.matrix(genotypes[!missidx,])
        covariates<-covariates[!missidx,]
        n<-length(phenotype)
	Z<-genotypes
	nZ<-ncol(Z)
	MAF<-colMeans(Z, na.rm = TRUE)/2
	if (is.function(mainweights)) {
		weights1<-mainweights(MAF)
	} else {
                if(!is.vector(mainweights)) stop("Check the class of your variable mainweights: should be function or vector!")
                if(length(mainweights)!=nZ) stop("Number of variants inconsistent between genotypes and mainweights. Check your data...")
		weights1<-mainweights
	}
	if (is.function(interweights)) {
		weights2<-interweights(MAF)
	} else {
                if(!is.vector(interweights)) stop("Check the class of your variable interweights: should be function or vector!")
                if(length(interweights)!=nZ) stop("Number of variants inconsistent between genotypes and interweights. Check your data...")
		weights2<-interweights
	}
	# impute missing genotypes
	nZ<-ncol(Z)
	for(pp in 1:nZ) {
		IDX<-which(is.na(Z[,pp]))
		if(length(IDX)>0) {
			if(binomialimpute) { # random imputation based on binomial distribution
				Z[IDX,pp]<-rbinom(length(IDX), 2, MAF[pp])
			} else { # default: set missing to 0
				Z[IDX,pp]<-0
			}
		}
	}
	# now Z should not have any missing values
	G<-t(t(Z)*weights1)
	EG<-t(t(Z)*weights2)*(covariates[,1]-mean(covariates[,1]))

	# null model
	X<-cbind(rep(1,n),as.matrix(covariates))
	tmpdata<-data.frame(phenotype,id=rep(1,n),covariates,G)
	nc<-ncol(covariates)
	if(nc==1) {
		fixedmod<-"phenotype ~ covariates"
	} else {
		fixedmod<-paste("phenotype ~", paste(names(tmpdata)[3:(2+nc)],collapse=" + "))
	}
	if(nZ==1) {
		lmmod<-paste(fixedmod,"G",sep=" + ")
	} else {
		lmmod<-paste(fixedmod, paste(names(tmpdata)[(nc+3):(2+nc+nZ)],collapse=" + "), sep=" + ")
	}
	X2<-cbind(X,as.matrix(G))
	nullmod2<-glm(as.formula(lmmod),data=tmpdata,family=fam)
	res2<-residuals(nullmod2, type="response")
	sigma2<-ifelse(fam=="binomial", 1, var(residuals(nullmod2, "pearson"))*(n-1)/n)
	Q2<-sum((t(EG) %*% res2)^2)/sigma2
	binvar2<-nullmod2$family$var(as.numeric(nullmod2$fitted))
	EGbEG<-t(EG*binvar2) %*% EG
	X2bEG<-t(X2*binvar2) %*% EG
	X2bX2<-t(X2*binvar2) %*% X2
	eig2<-eigen(EGbEG - t(X2bEG) %*% ginv(X2bX2) %*% X2bEG, symmetric=T, only.values=T)
	evals2<-eig2$values[zapsmall(eig2$values)>0]
	pINT_FIX<-pchisqsum(Q2,rep(1,length(evals2)),evals2,lower.tail=F,method="sad")
	return(pINT_FIX)
}

INT_RAN<-function(phenotype, genotypes, covariates, mainweights=wuweights, interweights=wuweights, family="gaussian", binomialimpute=FALSE) {
	fam<-try(match.arg(family, c("gaussian", "binomial")))
	if(class(fam)=="try-error") {
		print(family)
		stop("The argument 'family' has to be either gaussian or binomial")
	}
	# phenotype is a vector, genotypes is a matrix, covariates is a matrix
        if(class(phenotype)!="numeric" && class(phenotype)!="integer") stop("phenotype should be a numeric vector!")
        n<-length(phenotype)
        if(is.data.frame(genotypes)) genotypes<-as.matrix(genotypes)
        if(!is.matrix(genotypes)) stop("genotypes should be a matrix!")
        if(nrow(genotypes)!=n) stop("Number of individuals inconsistent between phenotype and genotypes. Check your data...")
        covariates<-as.matrix(covariates)
        if(nrow(covariates)!=n) stop("Number of individuals inconsistent between phenotype and covariates. Check your data...")
        missidx<-is.na(phenotype) | apply(is.na(covariates), 1, any)
        phenotype<-phenotype[!missidx]
        genotypes<-as.matrix(genotypes[!missidx,])
        covariates<-covariates[!missidx,]
        n<-length(phenotype)
	Z<-genotypes
	nZ<-ncol(Z)
	MAF<-colMeans(Z, na.rm = TRUE)/2
	if (is.function(mainweights)) {
		weights1<-mainweights(MAF)
	} else {
                if(!is.vector(mainweights)) stop("Check the class of your variable mainweights: should be function or vector!")
                if(length(mainweights)!=nZ) stop("Number of variants inconsistent between genotypes and mainweights. Check your data...")
		weights1<-mainweights
	}
	if (is.function(interweights)) {
		weights2<-interweights(MAF)
	} else {
                if(!is.vector(interweights)) stop("Check the class of your variable interweights: should be function or vector!")
                if(length(interweights)!=nZ) stop("Number of variants inconsistent between genotypes and interweights. Check your data...")
		weights2<-interweights
	}
	# impute missing genotypes
	nZ<-ncol(Z)
	for(pp in 1:nZ) {
		IDX<-which(is.na(Z[,pp]))
		if(length(IDX)>0) {
			if(binomialimpute) { # random imputation based on binomial distribution
				Z[IDX,pp]<-rbinom(length(IDX), 2, MAF[pp])
			} else { # default: set missing to 0
				Z[IDX,pp]<-0
			}
		}
	}
	# now Z should not have any missing values
	G<-t(t(Z)*weights1)
	EG<-t(t(Z)*weights2)*(covariates[,1]-mean(covariates[,1]))

	# null model
	X<-cbind(rep(1,n),as.matrix(covariates))
	tmpdata<-data.frame(phenotype,id=rep(1,n),covariates,G)
	nc<-ncol(covariates)
	if(nc==1) {
		fixedmod<-"phenotype ~ covariates"
	} else {
		fixedmod<-paste("phenotype ~", paste(names(tmpdata)[3:(2+nc)],collapse=" + "))
	}
	if(nZ==1) {
		randommod<-"~ 0 + G"
	} else {
		randommod<-paste("~ 0 +", paste(names(tmpdata)[(nc+3):(2+nc+nZ)],collapse=" + "))
	}
	nullmod1<-glmmPQL(fixed=as.formula(fixedmod),data=tmpdata,random=list(id=pdIdent(as.formula(randommod))),family=fam,niter=300,verbose=F)
	phat<-as.numeric(predict(nullmod1, data=tmpdata, type="response", level=1))
	DELTA<-if(fam=="binomial") diag(1/phat/(1-phat)) else diag(n)
	res1<-tmpdata$phenotype - phat
	varc<-VarCorr(nullmod1)
	deltai<-if(fam=="binomial") phat*(1-phat)/as.numeric(varc[nrow(varc),1]) else rep(1/as.numeric(varc[nrow(varc),1]),n)
	Gdeltai<-G*deltai
	Di<-solve(t(G)%*%(Gdeltai)+diag(rep(1/as.numeric(varc[1,1]),nZ)))
	XSX<-t(X)%*%(X*deltai)-t(X)%*%Gdeltai%*%Di%*%(t(Gdeltai)%*%X)
	XSEG<-t(X)%*%(EG*deltai)-t(X)%*%Gdeltai%*%Di%*%(t(Gdeltai)%*%EG)
	EGSEG<-t(EG)%*%(EG*deltai)-t(EG)%*%Gdeltai%*%Di%*%(t(Gdeltai)%*%EG)
	Q1<-sum((t(EG) %*% res1 / as.numeric(varc[nrow(varc),1]))^2)
	eig1<-eigen(EGSEG-t(XSEG)%*%solve(XSX,XSEG), symmetric=T, only.values=T)
	evals1<-eig1$values[zapsmall(eig1$values)>0]
	pINT_RAN<-pchisqsum(Q1,rep(1,length(evals1)),evals1,lower.tail=F,method="sad")
	return(pINT_RAN)
}

JOINT<-function(phenotype, genotypes, covariates, mainweights=wuweights, interweights=wuweights, family="gaussian", binomialimpute=FALSE, rho=seq(0, 1, by=0.1), B=10000) {
	fam<-try(match.arg(family, c("gaussian", "binomial")))
	if(class(fam)=="try-error") {
		print(family)
		stop("The argument 'family' has to be either gaussian or binomial")
	}
	# phenotype is a vector, genotypes is a matrix, covariates is a matrix
        if(class(phenotype)!="numeric" && class(phenotype)!="integer") stop("phenotype should be a numeric vector!")
        n<-length(phenotype)
        if(is.data.frame(genotypes)) genotypes<-as.matrix(genotypes)
        if(!is.matrix(genotypes)) stop("genotypes should be a matrix!")
        if(nrow(genotypes)!=n) stop("Number of individuals inconsistent between phenotype and genotypes. Check your data...")
        covariates<-as.matrix(covariates)
        if(nrow(covariates)!=n) stop("Number of individuals inconsistent between phenotype and covariates. Check your data...")
        missidx<-is.na(phenotype) | apply(is.na(covariates), 1, any)
        phenotype<-phenotype[!missidx]
        genotypes<-as.matrix(genotypes[!missidx,])
        covariates<-covariates[!missidx,]
        n<-length(phenotype)
	Z<-genotypes
	nZ<-ncol(Z)
	MAF<-colMeans(Z, na.rm = TRUE)/2
	if (is.function(mainweights)) {
		weights1<-mainweights(MAF)
	} else {
                if(!is.vector(mainweights)) stop("Check the class of your variable mainweights: should be function or vector!")
                if(length(mainweights)!=nZ) stop("Number of variants inconsistent between genotypes and mainweights. Check your data...")
		weights1<-mainweights
	}
	if (is.function(interweights)) {
		weights2<-interweights(MAF)
	} else {
                if(!is.vector(interweights)) stop("Check the class of your variable interweights: should be function or vector!")
                if(length(interweights)!=nZ) stop("Number of variants inconsistent between genotypes and interweights. Check your data...")
		weights2<-interweights
	}
	# impute missing genotypes
	nZ<-ncol(Z)
	for(pp in 1:nZ) {
		IDX<-which(is.na(Z[,pp]))
		if(length(IDX)>0) {
			if(binomialimpute) { # random imputation based on binomial distribution
				Z[IDX,pp]<-rbinom(length(IDX), 2, MAF[pp])
			} else { # default: set missing to 0
				Z[IDX,pp]<-0
			}
		}
	}
	# now Z should not have any missing values
	G<-t(t(Z)*weights1)
	EG<-t(t(Z)*weights2)*(covariates[,1]-mean(covariates[,1]))

	# null model
	X<-cbind(rep(1,n),as.matrix(covariates))
	tmpdata<-data.frame(phenotype,id=rep(1,n),covariates,G)
	nc<-ncol(covariates)
	if(nc==1) {
		fixedmod<-"phenotype ~ covariates"
	} else {
		fixedmod<-paste("phenotype ~", paste(names(tmpdata)[3:(2+nc)],collapse=" + "))
	}
	pJOINT<-pJOINTmin<-pJOINTrho<-pJOINTps<-pJOINTinfo<-NULL
	nullmod3<-glm(as.formula(fixedmod),data=tmpdata,family=fam)
	sigma2_3<-ifelse(fam=="binomial", 1, var(residuals(nullmod3, "pearson"))*(n-1)/n)
	res3<-residuals(nullmod3, type="response")/sqrt(sigma2_3)
	binvar3<-nullmod3$family$var(as.numeric(nullmod3$fitted))
	VV<-diag(binvar3)-(X*binvar3) %*% ginv(t(X*binvar3) %*% X) %*% t(X*binvar3)
	JOINT<-rareGE_joint_mc(res3, VV, G, EG, rho=rho, B=B)
	pJOINT<-JOINT$actualp
	pJOINTmin<-JOINT$minp
	pJOINTrho<-JOINT$rho
	pJOINTps<-JOINT$ps
	pJOINTinfo<-JOINT$info
	return(list(pJOINT=pJOINT,pJOINTmin=pJOINTmin,pJOINTrho=pJOINTrho,pJOINTps=pJOINTps,pJOINTinfo=pJOINTinfo))
}
