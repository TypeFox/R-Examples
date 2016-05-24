powerGG <- function(n, power, model, caco=0.5, alpha=0.05, alpha1)
{
        power.interaction.error <- function(xx)
        {
                warning(xx)
                return(1)
        }
        ww <- 0
        if(missing(n)+missing(power)!=1)stop("need to specify exactly one of n and power")
        if(!missing(n))if(n<10)ww <- power.interaction.error("combined sample size n is too small; should be at least 10")
        if(!missing(power))if(power<=0 || power>=1)ww <- power.interaction.error("required power should be between 0 and 1")
        if(model$pGene1<=0 || model$pGene1>=1)ww <- power.interaction.error("gene 1 probability pGene1 should be between 0 and 1")
        if(model$pGene2<=0 || model$pGene2>=1)ww <- power.interaction.error("gene 2 probability pGene1 should be between 0 and 1")
        if(length(model$beta.LOR)!=3)ww <- power.interaction.error("Need to specify three log odds ratios: [1] gene1 [2] gene2 [3] GxG")
        if(caco<=0 || caco>=1)ww <- power.interaction.error("case fraction caco should be between 0 and 1")
        if(model$nSNP<3)ww <- power.interaction.error("nSNP - number of SNPs tested - needs to be at least 3")
        if(alpha<=0 || alpha>=1)ww <- power.interaction.error("type 1 error alpha should be between 0 and 1")
        if(alpha1<=0 || alpha1>1)ww <- power.interaction.error("phase 1 alpha (alpha1) should be larger than 0 and at most 1")
        if(ww==1)stop("input errors")
        if(model$prev<=0 || model$prev>=1){
                warning(paste("disease prevalence p needs to be between 0 and 1, you specified",model$prev,"\nreset to 0.01: rare disease"))
                model$prev <- 0.01
        }
# if targ <1 the target is a power, if targ>0 the target is a samplesize
        if(missing(n)){
                targ <- 0
        } else {
                targ <- 1
        }

##############################################################
	doonepower <- function (b, maf, n, caco, nsnps, alpha1, crit)
	{
		cc <- round(caco*n)
		cc <- c(cc,n-cc)
		zpower <- function (n1,n2,p1,p2, alpha)
		{
			# this routine computes the power of finding a difference
			# between two binomial samples - sample i (i=1,2) with sample
			# size ni and success pi. The Type 1 error (2 sided) is alpha.
			# it's a slight generalization of the R function
			# power.prop.test() for different sample sizes.
			zalpha <- qnorm(alpha/2, lower.tail=FALSE)
			q1 <- 1-p1
			q2 <- 1-p2
			pbar <- (p1+p2)/2
			qbar <- (q1+q2)/2
			r <- n2/n1
			top <- sqrt(r*(p1-p2)*(p1-p2)*n1)-zalpha*sqrt((r+1)*pbar*qbar)
			bottom <- sqrt(r*p1*q1+p2*q2)
			pnorm(top/bottom)
		}
		# the design matrix for the two different SNPs and their interaction
		desmat <- matrix(c(1,1,1,1,0,1,0,1,0,0,1,1,0,0,0,1),ncol=4)

		# the fitted values, exponentiated
		p00 <- exp(desmat%*%as.matrix(b))

		# P(Y=0|X,b)
		q00 <- 1/(1+p00)

		# P(Y=1|X,b)
		p00 <- 1 - q00

		# now convert to conditional probabilities for case-control sampling
		q00 <- q00 * (1-maf[1])^(1-desmat[,2])*maf[1]^desmat[,2]
		p00 <- p00 * (1-maf[1])^(1-desmat[,2])*maf[1]^desmat[,2]
		q00 <- q00 * (1-maf[2])^(1-desmat[,3])*maf[2]^desmat[,3]
		p00 <- p00 * (1-maf[2])^(1-desmat[,3])*maf[2]^desmat[,3]
		q00 <- q00/sum(q00)
		p00 <- p00/sum(p00)

		# the power that SNP1 will be selected at stage 1
		z1 <- zpower(cc[1], cc[2], p00[2] + p00[4], q00[2] + q00[4], alpha1)

		# the power that SNP2 will be selected at stage 1
		z2 <- zpower(cc[1], cc[2], p00[3] + p00[4], q00[3] + q00[4], alpha1)

		# the power that both SNPs will be selected at stage 1
		zz <- z1*z2

		# expected standard error for the coefficient of interaction
		v <- sqrt((sum(1/q00)/cc[2]+sum(1/p00)/cc[1]))

		# v-caseonly
		vcaseonly <- sqrt(sum(1/p00)/cc[1])

		# normalized coefficient of the interaction parameter
		v <- abs(b[4]/v)
		vcaseonly <- abs(b[4]/vcaseonly)

		# the number of other (not interesting) SNPs
		nother <- nsnps-2

		# the plausible range of the number of other SNPs being selected
		u1 <- u3 <- 0
		qb <- (qbinom(.0000001,nother,alpha1)-10):(qbinom(1-.0000001,nother,alpha1)+10)
		qb <- qb[qb>=0]

		# the following lines yield vectors with different values for different qb
		# the probabilities of seeing that many other SNPs
		db <- dbinom(qb,nother,alpha1)

		# tails - just in case the probabilities don
		ub <- 1-sum(db)
		db[c(1,length(db))] <- db[c(1,length(db))]+ub/2
 
		# number of interactions involved in the selected SNPs
		nu <- (qb+2)*(qb+1)/2

		# Bonferonni alpha for each interaction
		alpha2 <- -qnorm(0.5*crit/nu)
    
		# power of identifying the "correct" interaction, if SNP1 and SNP2 survived stage 1
		p1 <- pnorm(-alpha2 - v) + 1 - pnorm(alpha2 - v)
		p1caseonly <- pnorm(-alpha2 - vcaseonly) + 1 - pnorm(alpha2 - vcaseonly)
  
		# integrating out: zz=P(snps survived stage 1), 
		# p1=P(identified|snps survived stage 1 and qb identified)
		# db=p(qb identified)
		u1 <- sum(zz*p1*db)
		u3 <- sum(zz*p1caseonly*db)
		c(u1,u3)
	}
#############################################################
	maf <- c(model$pGene1,model$pGene2)
	crit <- alpha
	nsnps <- c(model$nSNP)
	getbb <- function(maf,prev,bb)
	{
		getprev <- function(x,bb,maf,prev)
		{
			bb <- c(x,bb)
			pp <- bb[1]+bb[2]*c(0,1,0,1)+bb[3]*c(0,1,1,1)+c(0,0,0,bb[4])
			pp <- exp(pp)/(1+exp(pp))
			dd <- c((1-maf[1])*(1-maf[2]),maf[1]*(1-maf[2]),
				maf[2]*(1-maf[1]),maf[1]*maf[2])
			(sum(dd*pp)-prev)^2
		}	
		bb1 <-	nlminb(0,getprev,bb=bb,maf=maf,prev=prev)$par
		c(bb1,bb)	
	}
	bb <- getbb(maf,model$prev,model$beta.LOR)
	if(targ==1) {
		 po <- doonepower(bb, maf, n, caco, nsnps, alpha1, crit)
	}
	else{
		po <- c(NA,NA)
		for(j in 1:2){
			nmin <- 10
			pmin <- doonepower(bb,maf,nmin,caco,nsnps,alpha1,crit)[j]
			nmax <- 10^10
			pmax <- doonepower(bb,maf,nmax,caco,nsnps,alpha1,crit)[j]
			if(pmin[1]>=power){
				po[j] <- nmin
			}
			if(pmax[1]<power){
				po[j] <- -1
			}
			if(is.na(po[j]))for(ii in 1:30)if(nmax-nmin>1){
				nmid <- floor((nmax+nmin)/2)
				pmid <- doonepower(bb,maf,nmid,caco,nsnps,alpha1,crit)[j]
				if(pmid[1]>=power){
					nmax <- nmid
					pmax <- pmid
				} else {
					nmin <- nmid
					pmin <- pmid
				}
			}
			if(pmin[1]>=power){
				po[j] <- nmin
			}
			else{
				po[j] <- nmax
			}
		}
	}
	out <- data.frame(matrix(po,ncol=2,nrow=1))
	names(out) <- c("case.control","case.only")
	if(targ==1)row.names(out) <- "power"
	else row.names(out) <- "sample.size"
	out
}
powerGE <- function(n, power, model, caco=0.5, alpha=0.05, alpha1, maintain.alpha = TRUE)
{
	power.interaction.error <- function(xx)
	{
		warning(xx)
		return(1)	
	}
	ww <- 0
	if(missing(n)+missing(power)!=1)stop("need to specify exactly one of n and power")
	if(!missing(n))if(n<10)ww <- power.interaction.error("combined sample size n is too small; should be at least 10")
	if(!missing(power))if(power<=0 || power>=1)ww <- power.interaction.error("required power should be between 0 and 1")
	if(model$pGene<=0 || model$pGene>=1)ww <- power.interaction.error("gene probability pGene should be between 0 and 1")
	if(model$pEnv<=0 || model$pEnv>=1)ww <- power.interaction.error("binary enviornment  probability pGene should be between 0 and 1")
	if(length(model$beta.LOR)!=3)ww <- power.interaction.error("Need to specify three log odds ratios: [1] gene [2] enviornment [3] GxE")
	if(caco<=0 || caco>=1)ww <- power.interaction.error("case fraction caco should be between 0 and 1")
	if(model$nSNP<1)ww <- power.interaction.error("nSNP - number of SNPs tested - needs to be at least 1")
	if(model$orGE<=0)ww <- power.interaction.error("orGE - odds ratio between Gene and Environment needs to be larger than 0")
	if(alpha<=0 || alpha>=1)ww <- power.interaction.error("type 1 error alpha should be between 0 and 1")
	if(alpha1<=0 || alpha1>1)ww <- power.interaction.error("phase 1 alpha (alpha1) should be larger than 0 and at most 1")
	if(ww==1)stop("input errors")
	if(model$prev<=0 || model$prev>=1){
		warning(paste("disease prevalence p needs to be between 0 and 1, you specified",model$prev,"\nreset to 0.01: rare disease"))
		model$prev <- 0.01
	}
# if targ <1 the target is a power, if targ>0 the target is a samplesize
	if(missing(n)){
		targ <- 0
	} else {
		targ <- 1
	}

	prev <- model$prev
	pGene <- model$pGene
	pEnv <- model$pEnv
	beta.LOR <- model$beta.LOR
	orGE <- model$orGE
	nSNP <- model$nSNP
	sampleall <- matrix(NA,nrow=5,ncol=3)
	powerall <- matrix(NA,nrow=5,ncol=3)
	alphaall <- matrix(NA,nrow=5,ncol=3)
	# library(mvtnorm)
	# library(pwr)

	maintain <- 0
	if(maintain.alpha){
		if(orGE==1) maintain <- 1
		else maintain <- 2
	}
     
##################################################################################################
###  Main Program
### Power Comparisons of Two Stage Approaches 
### Input  targ: if targ>10 this is a sample size, if 0<targ<1 this is a power   
###           n:  sample sizes
### Input: model
###              prev: Disease prevalence 
###             pGene: P(G=1)
###              pEnv: P(E=1)
###          beta.LOR: log-odds ratio of the main effects of G and E, and the interaction effect of GxE. 
###              orGE: odds ratio of G and E in the population 
###              nSNP: # of SNPs tested
###        caco: fraction of cases and controls
###       alpha: Overall type I error (default=0.05)
###      alpha1: The alpha level at the screening stage. 
### Output: List two elements
###        [[1]]: 5x3 matrix with row (no screen, marg, corr, cocktail and chisquare) 
###               and column (cc, co, EB)   power or sample sizes
###        [[2]]: stage I alpha level that matches expected P-value.  If it is specified, otherwise it will be the same as the specified alpha1
###               if targ is n this one is 5x1, if targ is power this one is 5x3
###################################################################################################

### supplemental functions ###
      
###################################################################################################
	cond.p <- function(beta,G,E){
	####
	## calculate the logit disease probability given G,E and GxE
	## exp(beta*(G,E,GxE))/(1+exp(beta*(G,E,GxE)))
	####
		Z <- c(1,G,E,G*E)
		temp <- exp(sum(Z*beta))
		logit.p <- temp/(1+temp)
		return(logit.p)
	}
###################################################################################################
	cond.p2 <- function(beta,caco,prev,G,E){
	#####
	### Disease "probability" in the case-control sample
	### There is a shift in intercept due to sampling
	### caco: proportion of cases in the sample
	### prev: disease prevalence
	####
		Z <- c(1,G,E,G*E)
		beta1 <- beta
		beta1[1] <- log(caco/(1-caco)) - log(prev/(1-prev)) + beta[1]
		temp <- exp(sum(Z*beta1))
		logit.p <- temp/(1+temp)
		return(logit.p)
	}
###################################################################################################
	prob.GE0 <- function(a,pGene,pEnv,orGE){
	#####
	### This is the function to calculate a = Pr(G=1,E=1) given 
	### OR(G,E)=orGE and marginal probability of G and E
	## pGene: allele frequency of G
	## pEnv: frequency of exposure
	## orGE: Odds ratio of G and E; if orGE=1, G and E are independent
	## a: Pr(G=0, E=0)
	###
		obj <- (orGE-1)*a^2 - (2*orGE-1 - (orGE-1)*(pEnv+pGene))*a + orGE*(pEnv-1)*(pGene-1)
		return(obj^2)
	}
###################################################################################################
	marg.fcn <- function(gamma,caco,cov,cov1){
	####
	### Marginal association analysis
	###
		prob.g <- rowSums(cov)
		prob.g1 <- rowSums(cov1)
		temp1 <- temp2 <- 0
		for (i in 1:2){
			temp <- exp(sum(gamma*c(1,(i-1))))
			temp <- temp/(1+temp)
			temp1 <- temp1 + temp*prob.g[i]
			temp2 <- temp2 + (i-1)*temp*prob.g[i]
		}
		eqn1 <- caco - temp1
		eqn2 <- caco*prob.g1[2] - temp2
		eqn <- eqn1^2 + eqn2^2
		return(eqn)
	}
###################################################################################################
	marg.var <- function(gamma,cov){
	####
	### Calculation of variance for marginal association
	####
		dis.G <- vector()
		for (i in 1:2){
		      temp <- exp(sum(gamma*c(1,(i-1))))
		      dis.G[i] <- temp/(1+temp)
		}
	### Information matrix
		var0 <- matrix(0,2,2)
		prob.g <- rowSums(cov)
		for (i in 1:2){
		      tmp <- dis.G[i]*(1-dis.G[i])*prob.g[i]
		      var0[1,1] <- var0[1,1] + tmp
		      var0[1,2] <- var0[1,2] + (i-1)*tmp
		      var0[2,2] <- var0[2,2] + (i-1)^2*tmp
		}
		var0[2,1] <- var0[1,2]
	### Hessian matrix
		var1 <- solve(var0)
		return(var1)
	}
###################################################################################################
	corr.fcn <- function(gamma,cov){
	###
	### calculation of correlation
	### 
		prob.g <- rowSums(cov)
		prob.e <- colSums(cov)
		temp1 <- temp2 <- 0
		for (i in 1:2){
			temp <- exp(sum(gamma*c(1,(i-1))))
			temp <- temp/(1+temp)
			temp1 <- temp1 + temp*prob.e[i]
			temp2 <- temp2 + (i-1)*temp*prob.e[i]
		}
		eqn1 <- prob.g[2] - temp1
		eqn2 <- cov[2,2] - temp2
		eqn <- eqn1^2 + eqn2^2
		return(eqn)
	}
###################################################################################################
	corr.var <- function(gamma,cov){
	###      
	### calculation of variance for correlation screening
	### 
		dis.G <- vector()
		for (i in 1:2){
			temp <- exp(sum(gamma*c(1,(i-1))))
			dis.G[i] <- temp/(1+temp)
		}
	### Information matrix
		prob.e <- apply(cov,2,sum)
		var0 <- matrix(0,2,2)
		for (i in 1:2){
			tmp <- dis.G[i]*(1-dis.G[i])*prob.e[i]
			var0[1,1] <- var0[1,1] + tmp
			var0[1,2] <- var0[1,2] + (i-1)*tmp
		}
		var0[2,1] <- var0[1,2]
		var0[2,2] <- var0[1,2]
	### Hessian matrix
		var1 <- solve(var0)
		return(var1)
	}
###################################################################################################
	marg.corr.var <- function(gamma.mg,gamma.corr,cov,cov1,cov0,caco){
	###
	### Calculate the covariance between marginal and correlation screening
	###
		Hess <- matrix(0,4,4)
		Hess[1:2,1:2] <- marg.var(gamma.mg,cov)
		Hess[3:4,3:4] <- corr.var(gamma.corr,cov)
	#### The middle part of the sandwich variance estimator
		score <- function(D, G, E){
			l <- vector()
			P.G <- 1/(1+exp(-sum(c(1,G)*gamma.mg)))
			P.E <- 1/(1+exp(-sum(c(1,E)*gamma.corr)))
			l[1] <- (D-P.G)
			l[2] <- G*(D-P.G)
			l[3] <- (G-P.E)
			l[4] <- E*(G-P.E)
			return(l)
		}
	###
		P.D<-c((1-caco), caco)
		info<-matrix(0,4,4)
		for (k in 1:2) for (i in 1:2) for (j in 1:2){
			if (k==1) {
				cc <- cov0[i,j]
			} else {
				cc <- cov1[i,j]
			}
			tscore <- score((k-1),(i-1), (j-1))
			info <- info + (tscore %*% t(tscore))*P.D[k]*cc
		}
		var.all <- Hess %*% info %*% Hess
		return(var.all)
	}
###################################################################################################
	cc.var <- function(beta,cov){
		temp <- rep(0,4)
		for (i in 1:2) for (j in 1:2){
			dis.GE <- cond.p(beta,(i-1), (j-1)) * (1-cond.p(beta,(i-1), (j-1)))
			tmp <- dis.GE*cov[i,j]
			temp[1] <- temp[1] + tmp
			temp[2] <- temp[2] + (i-1)*tmp  ## E{G}
			temp[3] <- temp[3] + (j-1)*tmp  ## E{E}
			temp[4] <- temp[4] + (i-1)*(j-1)*tmp  ## E{GE}
		}
		var.cc <- matrix(0,4,4)
		var.cc[1,] <- temp
		var.cc[2,] <- c(temp[2],temp[2],temp[4],temp[4])
		var.cc[3,] <- c(temp[3],temp[4],temp[3], temp[4])
		var.cc[4,] <- rep(temp[4],4)
		var1 <- solve(var.cc)
		return(var1)
	}
###################################################################################################
	find.beta0 <- function(beta0,beta.LOR,prev,pGE){
	##############################################################################
	### This function is to calculate the intercept given
	### disease prevalence (prev), beta.OR(main effects and interaction effect, pGene, pEnv)
	################################################################################
		beta <- c(beta0,beta.LOR)
		prob0 <- 0
		for (i in 1:2) for (j in 1:2){
			prob0 <- prob0 + cond.p(beta,(i-1),(j-1))*pGE[i,j]
		}
		obj <- (prob0 - prev)^2
		return(obj)      
	}
###################################################################################################
### computations independent of power and n
###
### Joint distribution of G and E 
### pGE[1,1]=Pr(G=0,E=0); pGE[1,2]=Pr(G=0, E=1); pGE[2,1]=Pr(G=1,E=0)
### pGE[2,2] = Pr(G=1, E=1)
	pGE <- matrix(0,2,2)
	if (orGE==1) {
		res <- (1-pGene)*(1-pEnv)
	} else {
		res <- nlminb(0,prob.GE0, pGene=pGene, pEnv=pEnv, orGE=orGE,lower=0)$par
	}

	pGE[1,1] <- res                      ## G=0 E=0
	pGE[1,2] <- 1 - pGene - res          ## G=0 E=1
	pGE[2,1] <- 1 - pEnv - res           ## G=1 E=0
	pGE[2,2] <- pEnv + pGene + res - 1   ## G=1 E=1

### only the preparation calculations that do not depend on n

### Obtaining intercept based on disease prevalence, pGE and beta.OR.
	beta0 = nlminb(0, find.beta0, beta.LOR=beta.LOR, prev=prev, pGE=pGE)$par

### vector of coefficients including intercep
	beta <-  c(beta0, beta.LOR)
      
## cov1=Pr(G,E|D=1), cov0=Pr(G,E|D=0)
## cov =Pr(G,E) in the case-control sample
	cov <- matrix(0,2,2);cov1<-cov;cov0<-cov
	for (i in 1:2) for (j in 1:2){
		cov1[i,j] <- cond.p(beta,(i-1),(j-1))*pGE[i,j]/prev
		cov0[i,j] <- (1-cond.p(beta,(i-1),(j-1)))*pGE[i,j]/(1-prev)
		cov[i,j] <- cov1[i,j]*caco + cov0[i,j]*(1-caco)
	}
### Marginal Screening
	gamma.mg <- nlminb(c(0,0),marg.fcn,caco=caco,cov=cov,cov1=cov1)$par
	var.mgn <- marg.var(gamma.mg, cov)[2,2]
	z.mgn <- gamma.mg[2]/sqrt(var.mgn)

### Correlation Screening
	gamma.corr <- nlminb(c(0,0),corr.fcn,cov=cov)$par
	var.corrn <- corr.var(gamma.corr,cov)[2,2]
	z.corrn <- gamma.corr[2]/sqrt(var.corrn)

### covariance between marginal association and Correlation
	var1 <- marg.corr.var(gamma.mg,gamma.corr,cov,cov1,cov0,caco)

### Case-control GxE interaction effect
	beta.cc <- beta
	beta.cc[1] <- beta[1] + log(caco/(1-caco)) - log(prev/(1-prev))
	var.ccn <- cc.var(beta.cc,cov)[4,4]

### Case-only GxE interaction effect
	beta.co <- nlminb(c(0,0),corr.fcn,cov=cov1)$par
	var.con <- corr.var(beta.co,cov1)[2,2]/caco

### EB-estimator for GxE interaction
### calculation of thetaGE = logOR(G|E) in controls
	thetaGE <- nlminb(c(0,0),corr.fcn,cov=cov0)$par
	varGEn <- corr.var(thetaGE,cov0)[2,2]/(1-caco)

###############################
### stuff that depends on n ###
###############################
	computeonen <- function(n,z.mgn,z.corrn,var1,beta.cc,var.ccn,beta.co,var.con,thetaGE,varGEn,caco)
	{
### Marginal Screening
		z.mg <- z.mgn*sqrt(n)
		pvalue.mg<- (1-pnorm(abs(z.mg)))*2      
### Correlation Screening
		z.corr <- z.corrn*sqrt(n)
		pvalue.corr <- (1-pnorm(abs(z.corr)))*2      
### covariance between marginal association and Correlation
		par.var <- var1[c(2,4),c(2,4)]/n
		tt <- cov2cor(par.var)
### Cocktail approach
		cc <- max(abs(z.mg),abs(z.corr))
		pvalue.cocktail <- (1-pmvnorm(lower=c(-cc,-cc),upper=c(cc,cc),corr=tt)[1])
### chisquare screening approach      
		par.val <- as.matrix(c(z.mg,z.corr))
		zvalue <- t(par.val) %*% solve(tt) %*% par.val
		pvalue.chi <- 1-pchisq(zvalue, df=2)
### Case-control GxE interaction effect
		var.cc <- var.ccn/n
		cc.GXE <- abs(beta.cc[4]/sqrt(var.cc*n))      ### Effect size for cc test of GxE
### Case-only GxE interaction effect
		var.co <- var.con/n
		co.GXE <- beta.co[2]/sqrt(var.co*n*caco)   ### Effect size for co test of GxE
### EB-estimator for GxE interaction
### calculation of thetaGE = logOR(G|E) in controls
		varGE <- varGEn/n
		wt  <- thetaGE[2]^2/(thetaGE[2]^2 + var.cc)
		beta.EB  <- beta.cc[4]*wt + beta.co[2]*(1-wt)
		var.EB  <- var.co + varGE*(thetaGE[2]^2*(thetaGE[2]^2+3*var.cc)/(var.cc+thetaGE[2]^2)^2)^2
		eb.GXE  <- beta.EB/sqrt(var.EB*n)
		co.prob <- 1-pnorm(z.corr - z.mg)
		allps <- c(pvalue.mg,pvalue.corr,pvalue.cocktail, pvalue.chi)
		return(list(allps=allps,cc.GXE=cc.GXE,co.GXE=co.GXE,eb.GXE=eb.GXE,co.prob=co.prob))
	}

#############################################################################
### for a given n computes the power for one method #########################
### nlist arguments come from computeonen ###################################
### i1 is prescreening method (1-none, 2=marg 3-corr, 4=cocktail 5=chi-sq ###
### i2 is second stage 1=cc, 2=co, 3=eb #####################################
#############################################################################
	pown <- function(nlist,n,nSNP,alpha1,alpha,i1,i2)
	{
		alpha.expected <- nlist$allps
		cc.GXE <- nlist$cc.GXE
		co.GXE <- nlist$co.GXE
		eb.GXE <- nlist$eb.GXE
		powervec  <- c(0,0,0)
		final.alpha <- c(NA,alpha.expected)[i1]
		z.expected <- qnorm(1-alpha.expected/2)
		z.fixed <- qnorm(1-alpha1/2)
###### Calculation of the probability of the screening statistics greater than the threshold
		prob.larger <- c(1,(1+pnorm(-z.fixed-z.expected)-pnorm(z.fixed - z.expected)))[i1]
		no.tests <- ceiling(nSNP*c(1,rep(alpha1,4)))[i1]
		if(no.tests == 0)powwow <- 0
		else{
			if(i1==4){
				powervec[1] <- pwr.norm.test(d=cc.GXE,n=n,sig.level=alpha/no.tests)$power
				powervec[2] <- pwr.norm.test(d=co.GXE,n=(n*caco),sig.level=alpha/no.tests)$power
				powervec[3] <- pwr.norm.test(d=eb.GXE, n = n, sig.level=alpha/no.tests)$power
				powervec[2:3] <- nlist$co.prob*powervec[2:3]+(1-nlist$co.prob)*powervec[1]
				powwow <- powervec[i2]
			} else{
				if(i2==1) powwow <- pwr.norm.test(d=cc.GXE,n=n,sig.level=alpha/no.tests)$power
				if(i2==2) powwow <- pwr.norm.test(d=co.GXE,n=(n*caco),sig.level=alpha/no.tests)$power
				if(i2==3) powwow <- pwr.norm.test(d=eb.GXE, n = n, sig.level=alpha/no.tests)$power
			}
			powwow <- powwow*prob.larger
		}
		return(c(powwow,final.alpha,prob.larger))
	}

################# compute samplesize for given power ##################
	if(targ==0){
		problarg <- matrix(NA,nrow=5,ncol=3)
		for(i1 in 1:5)for(i2 in 1:3){
			if((maintain>0 && i2>1 && (i1==2 || i1==5))||(maintain==2 && i2==2)){
				sampleall[i1,i2] <- -1
				alphaall[i1,i2] <- NA
			} else {
				nmin <- 10
				tmp <- computeonen(nmin,z.mgn,z.corrn,var1,beta.cc,var.ccn,beta.co,var.con,thetaGE,varGEn,caco)
				pmin <- pown(tmp,nmin,nSNP,alpha1,alpha,i1,i2)
				nmax <- 10^10
				tmp <- computeonen(nmax,z.mgn,z.corrn,var1,beta.cc,var.ccn,beta.co,var.con,thetaGE,varGEn,caco)
				pmax <- pown(tmp,nmax,nSNP,alpha1,alpha,i1,i2)
				if(pmin[1]>=power){
					sampleall[i1,i2] <- nmin
					alphaall[i1,i2] <- pmin[2]
					problarg[i1,i2] <- pmin[3]
				}
				if(pmax[1]<power){
					sampleall[i1,i2] <- -1
					alphaall[i1,i2] <- NA
					problarg[i1,i2] <- NA
				}
				if(is.na(sampleall[i1,i2]))for(ii in 1:30)if(nmax-nmin>1){
					nmid <- floor((nmax+nmin)/2)
					tmp <- computeonen(nmid,z.mgn,z.corrn,var1,beta.cc,var.ccn,beta.co,var.con,thetaGE,varGEn,caco)
					pmid <- pown(tmp,nmid,nSNP,alpha1,alpha,i1,i2)
					if(pmid[1]>=power){
						nmax <- nmid
						pmax <- pmid
					} else {
						nmin <- nmid
						pmin <- pmid
					}
				}
				if(pmin[1]>=power){
					sampleall[i1,i2] <- nmin
					alphaall[i1,i2] <- pmin[2]
					problarg[i1,i2] <- pmin[3]
				} else {
					sampleall[i1,i2] <- nmax
					alphaall[i1,i2] <- pmax[2]
					problarg[i1,i2] <- pmax[3]
				}
			}
		}
		sampleall[sampleall<0] <- NA
		out <- list(sampleall,alphaall,problarg)
	} else {
################# compute power for given samplesize ##################
		nlist <- computeonen(n,z.mgn,z.corrn,var1,beta.cc,var.ccn,beta.co,var.con,thetaGE,varGEn,caco)
		alpha.expected <- nlist$allps
		cc.GXE <- nlist$cc.GXE
		co.GXE <- nlist$co.GXE
		eb.GXE <- nlist$eb.GXE
		final.alpha <- matrix(c(NA,alpha.expected),ncol=5,nrow=3)
		z.expected <- qnorm(1-alpha.expected/2)
		z.fixed <- qnorm(1-alpha1/2)
###### Calculation of the probability of the screening statistics greater than the threshold
  		prob.larger <- c(1,(1+pnorm(-z.fixed-z.expected)-pnorm(z.fixed - z.expected)))
		no.tests <- ceiling(nSNP*c(1,rep(alpha1,4)))
		for (i in 1:5){
			if (no.tests[i] == 0){
				powerall[i,] = 0
			} else {
				powerall[i,1] <- pwr.norm.test(d=cc.GXE,n=n,sig.level=alpha/no.tests[i])$power
				powerall[i,2] <- pwr.norm.test(d=co.GXE,n=(n*caco),sig.level=alpha/no.tests[i])$power
				powerall[i,3] <- pwr.norm.test(d=eb.GXE, n = n, sig.level=alpha/no.tests[i])$power
### Special treatment of cocktail
				if(i==4){
					powerall[4,2:3] <- nlist$co.prob*powerall[4,2:3]+(1-nlist$co.prob)*powerall[4,1]
				}
			}
		}
                powerall <- powerall*prob.larger
		out <- list(powerall,final.alpha,matrix(prob.larger,ncol=5,nrow=3))
	}

## reformat output ###
	for(i in 1:3){
		out[[i]] <- data.frame(matrix(out[[i]],nrow=5))
		row.names(out[[i]]) <- c("st1.no.filter","st1.marginal.screen","st1.correlation.screen","st1.cocktail.screen","st1.chi.square.screen")
		names(out[[i]]) <- c("st2.case.control","st2.case.only","ste2.empirical.Bayes")
		if(maintain.alpha && model$orGE!=1) out[[i]][,2] <- NA
		if(maintain.alpha) out[[i]][c(3,5),c(2,3)] <- NA
	}

	if(model$orGE!=1 & maintain.alpha==FALSE){
		cat(paste("\n==> G and E are correlated - case-only estimaters do not maintain the Type 1 error\n",
			"The use of these estimators in that situation is not recommended.\n",
			"The Empirical Bayes estimators protect the Type 1 error better when the G - E correlation\n",
			"is small. When the G - E correlation is large the case-control estimator is the only one\n",
			"that maintains the Type 1 error.\n\n"))
	}
	if(maintain.alpha==FALSE){
		cat(paste("\n==> the case-only estimator and the Empirical Bayes estimator after correlation screening,\n",       
			"or chi-square screening does not maintain the Type 1 error",
			"The use of these estimators in that situation is not recommended.\n",
			"Use of these estimators after marginal screening or cocktail screening does maintain the correct Type 1 error.\n\n"))
	}
	if(model$prev>0.1)cat(paste("\n==> Your disease probability is larger than 0.1. It is questionable whether the case-only\n",
			"and Empirical Bayes estimators maintain the Type 1 error.\n\n"))

	if(missing(power)){
		return(list(power=out[[1]],expected.p=out[[2]],prob.select=out[[3]]))
	} else {
		return(list(samplesize=out[[1]],expected.p=out[[2]],prob.select=out[[3]]))
	}
}
powerGWASinteraction <- function()
{
	cat(paste("This function is depreciated. use\npowerGE for powercalculations for GxE interactions\n",
		"powerGG for powercalculations for GxG interactions\n"))
	invisible()
}
