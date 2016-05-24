#' Fast score test for association
#' 
#' Fast score test for association 
#' between a trait and genetic polymorphism
#' 
#' When formula contains covariates, the traits is analysed using GLM and later 
#' residuals used when score test is computed for each of the SNPs in analysis. 
#' Coefficients of regression are reported for the quantitative trait.
#' 
#' For binary traits, odds ratios (ORs) are reportted. When adjustemnt is 
#' performed, first, "response" residuals are estimated after adjustment for 
#' covariates and scaled to [0,1]. Reported effects are approximately equal 
#' to ORs expected in logistic regression model. 
#' 
#' With no adjustment for binary traits, 1 d.f., the test is equivalent to the 
#' Armitage test.
#' 
#' This is a valid function to analyse GWA data, including X chromosome. For X chromosome,
#' stratified analysis is performed (strata=sex).
#' 
#' @param formula Formula describing fixed effects to be used in analysis, e.g. 
#' y ~ a + b means that outcome (y) depends on two covariates, a and b. 
#' If no covariates used in analysis, skip the right-hand side of the 
#' equation.
#' @param data An object of \code{\link{gwaa.data-class}}
#' @param snpsubset ndex, character or logical vector with subset of SNPs to run analysis on. 
#' If missing, all SNPs from \code{data} are used for analysis.
#' @param idsubset ndex, character or logical vector with subset of IDs to run analysis on. 
#' If missing, all people from \code{data/cc} are used for analysis.
#' @param strata Stratification variable. If provieded, scores are computed within strata and 
#' then added up.
#' @param trait.type "gaussian" or "binomial" or "guess" (later option guesses trait type)
#' @param times If more than one, the number of replicas to be used in derivation of 
#' empirical genome-wide significance. See \code{\link{emp.qtscore}}, which
#' calls qtscore with times>1 for details
#' @param quiet do not print warning messages
#' @param bcast If the argument times > 1, progress is reported once in bcast replicas
#' @param clambda If inflation facot Lambda is estimated as lower then one, this parameter 
#' controls if the original P1df (clambda=TRUE) to be reported in Pc1df, 
#' or the original 1df statistics is to be multiplied onto this "deflation" 
#' factor (clambda=FALSE). If a numeric value is provided, it is used as a correction factor.
#' @param propPs proportion of non-corrected P-values used to estimate the inflation factor Lambda,
#' passed directly to the \code{\link{estlambda}}
#' @param details when FALSE, SNP and ID names are not reported in the returned object
#' (saves some memory). This is experimental and will be not mantained anymore 
#' as soon as we achieve better memory efficiency for storage of SNP and ID names
#' (currently default R character data type used)
#' 
#' @return Object of class \code{\link{scan.gwaa-class}}
#' 
#' @references 
#' Aulchenko YS, de Koning DJ, Haley C. Genomewide rapid association using mixed model 
#' and regression: a fast and simple method for genome-wide pedigree-based quantitative 
#' trait loci association analysis. Genetics. 2007 177(1):577-85.
#' 
#' Amin N, van Duijn CM, Aulchenko YS. A genomic background based method for 
#' association analysis in related individuals. PLoS ONE. 2007 Dec 5;2(12):e1274. 
#' 
#' @author Yurii Aulchenko
#' 
#' @seealso \code{\link{mlreg}},
#' \code{\link{mmscore}},
#' \code{\link{egscore}},
#' \code{\link{emp.qtscore}},
#' \code{\link{plot.scan.gwaa}},
#' \code{\link{scan.gwaa-class}}
#' 
#' @examples 
#' data(srdta)
#' #qtscore with stratification
#' a <- qtscore(qt3~sex,data=srdta)
#' plot(a)
#' b <- qtscore(qt3,strata=phdata(srdta)$sex,data=srdta)
#' add.plot(b,col="green",cex=2)
#' # qtscore with extra adjustment
#' a <- qtscore(qt3~sex+age,data=srdta)
#' a
#' plot(a)
#' # compare results of score and chi-square test for binary trait
#' a1 <- ccfast("bt",data=srdta,snps=c(1:100))
#' a2 <- qtscore(bt,data=srdta,snps=c(1:100),trait.type="binomial")
#' plot(a1,ylim=c(0,2))
#' add.plot(a2,col="red",cex=1.5)
#' # the good thing about score test is that we can do adjustment...
#' a2 <- qtscore(bt~age+sex,data=srdta,snps=c(1:100),trait.type="binomial")
#' points(a2[,"Position"],-log10(a2[,"P1df"]),col="green")
#' 
#' @keywords htest
#' 
"qtscore" <-
		function(formula,data,snpsubset,idsubset,strata,trait.type="gaussian",
				times=1,quiet=FALSE,bcast=10,clambda=TRUE,propPs=1.0,details=TRUE) {
	if (!is(data,"gwaa.data")) {
		stop("wrong data class: should be gwaa.data")
	}
	checkphengen(data)
	if (!missing(snpsubset)) data <- data[,snpsubset]
	if (!missing(idsubset)) data <- data[idsubset,]
	if (missing(strata)) {nstra=1; strata <- rep(0,data@gtdata@nids)}
	ttargs <- c("gaussian","binomial","guess")
	if (!(match(trait.type,ttargs,nomatch=0)>0)) {
		out <- paste("trait.type argument should be one of",ttargs,"\n")
		stop(out)
	}

	if ( is(try(formula,silent=TRUE),"try-error") ) { 
		formula <- phdata(data)[[as(match.call()[["formula"]],"character")]] 
	}
	
	if (trait.type=="guess") {
#		if (!missing(data)) attach(data@phdata,pos=2,warn.conflicts=FALSE)
		if (is(formula,"formula")) {
			mf <- model.frame(formula,data,na.action=na.omit,drop.unused.levels=TRUE)
			y <- model.response(mf)
			if (isbinomial(y)) trait.type <- "binomial" else trait.type <- "gaussian"
		} else if (is(formula,"numeric") || is(formula,"integer") || is(formula,"double")) {
			y <- formula
			if (isbinomial(y)) trait.type <- "binomial" else trait.type <- "gaussian"
		} else {
			stop("formula argument must be a formula or one of (numeric, integer, double)")
		}
		warning(paste("trait type is guessed as",trait.type))
#		if (!missing(data)) detach(data@phdata)
	}
	if (trait.type=="gaussian") fam <- gaussian()
	if (trait.type=="binomial") fam <- binomial()
	
#	if (!missing(data)) attach(data@phdata,pos=2,warn.conflicts=FALSE)
	if (is(formula,"formula")) {
		mf <- model.frame(formula,data,na.action=na.omit,drop.unused.levels=TRUE)
		y <- model.response(mf)
		test.type(y,trait.type)
		desmat <- model.matrix(formula,mf)
		lmf <- glm.fit(desmat,y,family=fam)
		mids <- rownames(data@phdata) %in% rownames(mf)
		if (trait.type=="binomial") {
			resid <- residuals(lmf,type="response")
			resid <- (resid-min(resid))/(max(resid)-min(resid))
		} else {
			resid <- lmf$resid
		}
	} else if (is(formula,"numeric") || is(formula,"integer") || is(formula,"double")) {
		y <- formula
		test.type(y,trait.type)
		mids <- (!is.na(y))
		y <- y[mids]
#		if (trait.type=="binomial") {
#			y <- glm(y~rep(1,length(y)),family=binomial)$resid
#		}
		resid <- y
	} else {
		stop("formula argument should be a formula or a numeric vector")
	}
#	if (!missing(data)) detach(data@phdata)
	if (length(strata)!=nids(data)) stop("Strata variable and the data do not match in length")
	if (any(is.na(strata))) stop("Strata variable contains NAs")
	if (any(strata!=0)) {
		olev <- levels(as.factor(strata))
		nstra <- length(olev)
		tstr <- strata
		for (i in 0:(nstra-1)) tstr <- replace(tstr,(strata==olev[i+1]),i)
		strata <- tstr
		rm(tstr)
	}
	nstra <- length(levels(as.factor(strata)))
	
	tmeas <- as.logical(mids)
	strata <- strata[tmeas]
	
	if (any(tmeas == FALSE)) {
		if (!quiet) warning(paste(sum(!tmeas),"observations deleted due to missingness"))
		data <- data[tmeas,]
	}
	
#	if (trait.type=="binomial" & !is(formula,"formula")) bin<-1 else bin <- 0
	if (trait.type=="binomial") bin<-1 else bin<-0
	lenn <- nsnps(data);
	###out <- list()
	if (times>1) {pb <- txtProgressBar(style = 3)}
	for (j in c(1:(times+1*(times>1)))) {
		if (j>1) resid <- sample(resid,replace=FALSE)
#		if (old) {
			chi2 <- .C("qtscore_glob",as.raw(data@gtdata@gtps),as.double(resid),as.integer(bin),
					as.integer(data@gtdata@nids),as.integer(data@gtdata@nsnps), as.integer(nstra), 
					as.integer(strata), chi2 = double(10*data@gtdata@nsnps), PACKAGE="GenABEL")$chi2
#		} else {
#			gtNrow <- dim(data)[1]
#			gtNcol <- dim(data)[2]
#			chi2 <- .Call("iteratorGA", 
#					data@gtdata@gtps, 
#					as.integer(gtNrow), as.integer(gtNcol),
#					as.character("qtscore_glob"),   # Function
#					as.character("R"),  			# Output type
#					as.integer(2),      			# MAR
#					as.integer(1),					# Steps
#					as.integer(5),      			# nr of extra arguments
#					as.double(resid),
#					as.integer(bin),
#					as.integer(data@gtdata@nids),
#					as.integer(nstra), 
#					as.integer(strata), 
#					package="GenABEL")
#		}
		if (any(chromosome(data)=="X")) {
			ogX <- gtdata(data[,chromosome(data)=="X"])
			sxstra <- strata; sxstra[male(ogX)==1] <- strata[male(ogX)==1]+nstra
#			if (old) {
				chi2X <- .C("qtscore_glob",as.raw(ogX@gtps),as.double(resid),as.integer(bin),
						as.integer(nids(ogX)),as.integer(nsnps(ogX)), as.integer(nstra*2), 
						as.integer(sxstra), chi2 = double(10*nsnps(ogX)), PACKAGE="GenABEL")$chi2
#			} else {
#				chi2X <- .Call("iteratorGA", 
#						ogX@gtps, 
#						as.integer(ogX@nsnps), 	# nCol
#						as.integer(ogX@nids), 	# nRow
#						as.character("qtscore_glob"), # Function
#						as.character("R"),  	# Output type
#						as.integer(1),      	# MAR
#						as.integer(1),			# Steps
#						as.integer(5),      	# nr of arguments
#						as.double(resid),
#						as.integer(bin),
#						as.integer(nids),
#						as.integer(nstra*2), 
#						as.integer(sxstra), 
#						package="GenABEL")
#			}
			revec <- (chromosome(data)=="X")
			revec <- rep(revec,6)
			chi2 <- replace(chi2,revec,chi2X)
			rm(ogX,chi2X,revec);gc(verbose=FALSE)
		}
		if (j == 1) {
			chi2.1df <- chi2[1:lenn];
			chi2.1df[abs(chi2.1df+999.99)<1.e-8] <- 0 #NA
			###out$chi2.1df <- chi2.1df
			chi2.2df <- chi2[(lenn+1):(2*lenn)];
			chi2.2df[abs(chi2.2df+999.99)<1.e-8] <- 0 #NA
			actdf <- chi2[(2*lenn+1):(3*lenn)];
			actdf[abs(actdf+999.99)<1.e-8] <- 1.e-16 #NA
#			out$actdf <- actdf
			###out$chi2.2df <- chi2.2df
#			z0 <- chi2[(7*lenn+1):(8*lenn)];
#			z0[abs(z0+999.99)<1.e-8] <- 0 #NA
##			out$z0 <- z0
#			z2 <- chi2[(8*lenn+1):(9*lenn)];
#			z2[abs(z2+999.99)<1.e-8] <- 0 #NA
##			out$z2 <- z2
#			rho <- chi2[(9*lenn+1):(10*lenn)];
#			rho[abs(rho+999.99)<1.e-8] <- 0 #NA
##			rho <- abs(rho)
##			out$rho <- rho
			
			lambda <- list()
			if (is.logical(clambda)) {
				if (lenn<10) {
					warning("no. observations < 10; Lambda set to 1")
					lambda$estimate <- 1.0
					lambda$se <- NA
				} else {
					if (lenn<100) warning("Number of observations < 100, Lambda estimate is unreliable")
					lambda <- estlambda(chi2.1df,plot=FALSE,proportion=propPs)
					if (lambda$estimate<1.0 && clambda==TRUE) {
						warning("Lambda estimated < 1, set to 1")
						lambda$estimate <- 1.0
						lambda$se <- NA
					}
				}
			} else {
				if (is.numeric(clambda)) {
					lambda$estimate <- clambda
					lambda$se <- NA
				} else if (is.list(clambda)) {
					if (any(is.na(match(c("estimate","se"),names(clambda)))))
						stop("when clambda is list, should contain estimate and se")
					lambda <- clambda
					lambda$se <- NA
				} else {
					stop("clambda should be logical, numeric, or list")
				}
			}
			chi2.c1df <- chi2.1df/lambda$estimate
			
#			if (is.logical(clambda)) {
#				lambda$iz0 <- estlambda(z0*z0,plot=FALSE,proportion=propPs)$estimate 
#				lambda$iz2 <- estlambda(z2*z2,plot=FALSE,proportion=propPs)$estimate
#				if (clambda && lambda$iz0<1.0) {warning("z0 lambda < 1, set to 1");lambda$iz0<-1.0}
#				if (clambda && lambda$iz2<1.0) {warning("z2 lambda < 1, set to 1");lambda$iz2<-1.0}
#				chi2.c2df <- (z0*z0/lambda$iz0 + z2*z2/lambda$iz2 - 2.*z0*z2*rho/(sqrt(lambda$iz0*lambda$iz2)))/(1.- rho*rho)
#			} else {
#				if (is.list(clambda) && !any(is.na(match(c("estimate","iz0","iz2"),names(clambda))))) {
#					chi2.c2df <- (z0*z0/lambda$iz0 + z2*z2/lambda$iz2 - 2.*z0*z2*rho/(sqrt(lambda$iz0*lambda$iz2)))/(1.- rho*rho)
#				} else {
#					lambda$iz0 <- 1.0
#					lambda$iz2 <- 1.0
#					chi2.c2df <- chi2.2df
#				}
#			}
			effB <- chi2[(3*lenn+1):(lenn*4)]
			effB[abs(effB+999.99)<1.e-8] <- NA
			effAB <- chi2[(4*lenn+1):(lenn*5)]
			effAB[abs(effAB+999.99)<1.e-8] <- NA
			effBB <- chi2[(5*lenn+1):(lenn*6)]
			effBB[abs(effBB+999.99)<1.e-8] <- NA
			if (times>1) {
				pr.1df <- rep(0,lenn)
				pr.2df <- rep(0,lenn)
				pr.c1df <- rep(0,lenn)
				pr.c2df <- rep(0,lenn)
			}
		} else {
			th1 <- max(chi2[1:lenn])
			pr.1df <- pr.1df + 1*(chi2.1df < th1)
			pr.2df <- pr.2df + 1*(chi2.2df < max(chi2[(lenn+1):(2*lenn)]))
			pr.c1df <- pr.c1df + 1*(chi2.c1df < th1)
#			pr.c2df <- pr.c2df + 1*(chi2.c2df < th1)
#			if (!quiet && ((j-1)/bcast == round((j-1)/bcast))) {
			##				cat("\b\b\b\b\b\b",round((100*(j-1)/times),digits=2),"%",sep="")
#				cat(" ",round((100*(j-1)/times),digits=2),"%",sep="")
#				flush.console()
#			}
			if (!quiet && ((j-1)/bcast == round((j-1)/bcast))) {
				setTxtProgressBar(pb, (j-1)/times)
			}
		}
	}
	if (times > bcast) {setTxtProgressBar(pb, 1.0);cat("\n")}
	
	if (times>1) {
		P1df <- pr.1df/times
		P1df <- replace(P1df,(P1df==0),1/(1+times))
		P2df <- pr.2df/times
		P2df <- replace(P2df,(P2df==0),1/(1+times))
		Pc1df <- pr.c1df/times
		Pc1df <- replace(Pc1df,(Pc1df==0),1/(1+times))
		#out$Pc2df <- pr.c2df/times
		#out$Pc2df <- replace(out$Pc2df,(out$Pc2df==0),1/(1+times))
	} else {
		P1df <- pchisq(chi2.1df,1,lower.tail=F)
		P2df <- pchisq(chi2.2df,actdf,lower.tail=F)
		Pc1df <- NULL
		#out$Pc1df <- pchisq(chi2.c1df,1,lower=F)
		#out$Pc2df <- pchisq(chi2.c2df,2,lower=F)
	}
	###out$lambda <- lambda
	###out$effB <- effB
	###out$effAB <- effAB
	###out$effBB <- effBB
	###out$N <- chi2[(6*lenn+1):(lenn*7)]
	###if (details) {
	###	out$snpnames <- data@gtdata@snpnames
	###	out$idnames <- data@gtdata@idnames
	###}
	###out$map <- data@gtdata@map
	###out$chromosome <- data@gtdata@chromosome
	###out$formula <- match.call()
	###out$family <- paste("score test for association with trait type",trait.type)
	#class(out) <- "scan.gwaa"
	if (is.null(Pc1df)) {
		results <- data.frame(N=chi2[(6*lenn+1):(lenn*7)],
				effB = effB, se_effB = abs(effB)/sqrt(chi2.1df), chi2.1df = chi2.1df, P1df = P1df, 
				effAB=effAB, effBB=effBB, chi2.2df = chi2.2df, P2df = P2df,
				stringsAsFactors = FALSE)
	} else {
		results <- data.frame(N=chi2[(6*lenn+1):(lenn*7)],
				effB = effB, se_effB = abs(effB)/sqrt(chi2.1df), chi2.1df = chi2.1df, P1df = P1df, 
				Pc1df = Pc1df, 
				effAB=effAB, effBB=effBB, chi2.2df = chi2.2df, P2df = P2df,
				stringsAsFactors = FALSE)
	}
	rownames(results) <- snpnames(data)
	out <- new("scan.gwaa",
			results=results,
			annotation = annotation(data), 
			lambda = lambda,
			idnames = idnames(data), 
			call = match.call(), 
			family = trait.type
	) 
	out
}

