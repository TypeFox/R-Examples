#' Fast score test for association, corrected with PC
#' 
#' Fast score test for association between a trait and genetic polymorphism,
#' adjusted for possible stratification by principal components.
#' 
#' The idea of this test is to use genomic kinship matrix to first, derive axes
#' of genetic variation (principal components), and, second, adjust both trait
#' and genotypes onto these axes. Note that the diagonal of the kinship matrix
#' should be replaced (default it is 0.5*(1+F), and for EIGENSTRAT one needs
#' variance).  These variances are porduced by \code{\link{hom}} function (see
#' example).
#' 
#' The traits is first analysed using LM and with covariates as specified with
#' formula and also with axes of variation as predictors. Corrected genotypes
#' are defined as residuals from regression of genotypes onto axes (which are
#' orthogonal). Correlaton between corrected genotypes and phenotype is
#' computed, and test statistics is defined as square of this correlation times
#' (N - K - 1), where N is number of genotyped subjects and K is the number of
#' axes.
#' 
#' This test is defined only for 1 d.f.
#' 
#' @param formula Formula describing fixed effects to be used in analysis, e.g.
#' y ~ a + b means that outcome (y) depends on two covariates, a and b.  If no
#' covariates used in analysis, skip the right-hand side of the equation.
#' @param data An object of \code{\link{gwaa.data-class}}
#' @param snpsubset Index, character or logical vector with subset of SNPs to
#' run analysis on.  If missing, all SNPs from \code{data} are used for
#' analysis.
#' @param idsubset Index, character or logical vector with subset of IDs to run
#' analysis on.  If missing, all people from \code{data/cc} are used for
#' analysis.
#' @param kinship.matrix kinship matrix, as returned by \code{\link{ibs}}, Use
#' weight="freq" with \code{\link{ibs}} and do not forget to repalce the
#' diagonal with Var returned by \code{\link{hom}}, as shown in example!
#' @param naxes Number of axes of variation to be used in adjustment (should be
#' much smaller than number of subjects)
#' @param strata Stratification variable. If provieded, scores are computed
#' within strata and then added up.
#' @param times If more then one, the number of replicas to be used in
#' derivation of empirical genome-wide significance.
#' @param quiet do not print warning messages
#' @param bcast If the argument times > 1, progress is reported once in bcast
#' replicas
#' @param clambda If inflation facot Lambda is estimated as lower then one,
#' this parameter controls if the original P1df (clambda=TRUE) to be reported
#' in Pc1df, or the original 1df statistics is to be multiplied onto this
#' "deflation" factor (clambda=FALSE).  If a numeric value is provided, it is
#' used as a correction factor.
#' @param propPs proportion of non-corrected P-values used to estimate the
#' inflation factor Lambda, passed directly to the \code{\link{estlambda}}
#' @return Object of class \code{\link{scan.gwaa-class}}
#' @author Yurii Aulchenko
#' @seealso \code{\link{qtscore}}, \code{\link{mmscore}}, \code{\link{ibs}},
#' \code{\link{scan.gwaa-class}}
#' @references Price A. L. et al, Principal components analysis corrects for
#' stratification in genome-wide association studies. Nat Genet 38: 904-909.
#' @keywords htest
#' @examples
#' 
#' data(ge03d2c)
#' #egscore with stratification
#' gkin <- ibs(ge03d2c[,autosomal(ge03d2c)],w="freq")
#' #replace the diagonal with right elements
#' diag(gkin) <- hom(ge03d2c[,autosomal(ge03d2c)])$Var
#' a <- egscore(dm2~sex+age,data=ge03d2c,kin=gkin)
#' plot(a,df="Pc1df")
#' 
"egscore" <-
function(formula,data,snpsubset,idsubset,kinship.matrix,naxes=3,strata,
		times=1,quiet=FALSE,bcast=10,clambda=TRUE,propPs=1.0) {
  	if (!is(data,"gwaa.data")) {
		stop("wrong data class: should be gwaa.data")
  	}
	checkphengen(data)
	if (!missing(snpsubset)) data <- data[,snpsubset]
	if (!missing(idsubset)) data <- data[idsubset,]
	
	if (missing(strata)) {nstra=1; strata <- rep(0,data@gtdata@nids)}

	if (length(strata)!=data@gtdata@nids) stop("Strata variable and the data do not match in length")
	if (any(is.na(strata))) stop("Strata variable contains NAs")

	if ( is(try(formula,silent=TRUE),"try-error") ) { 
		formula <- phdata(data)[[as(match.call()[["formula"]],"character")]] 
	}
	
#	if (!missing(data)) attach(data@phdata,pos=2,warn.conflicts=FALSE)
	if (is(formula,"formula")) {
		mf <- model.frame(formula,data,na.action=na.omit,drop.unused.levels=TRUE)
		y <- model.response(mf)
		desmat <- model.matrix(formula,mf)
		lmf <- glm.fit(desmat,y)
		mids <- rownames(data@phdata) %in% rownames(mf)
		resid <- lmf$resid
	} else if (is(formula,"numeric") || is(formula,"integer") || is(formula,"double")) {
		y <- formula
		mids <- (!is.na(y))
		y <- y[mids]
		resid <- y
		if (length(unique(resid))==1) stop("trait is monomorphic")
	} else {
		stop("formula argument must be a formula or one of (numeric, integer, double)")
	}
#	if (!missing(data)) detach(data@phdata)
	if (length(strata)!=data@gtdata@nids) stop("Strata variable and the data do not match in length")
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
		if (!quiet) warning(paste(sum(!tmeas),"people (out of",length(tmeas),") excluded because they have trait or covariate missing\n"),immediate. = TRUE)
		data <- data[tmeas,]
	}

	tmp <- dim(kinship.matrix)
	if (tmp[1] == tmp[2]) {
		tmp <- t(kinship.matrix)
		diag(kinship.matrix) <- attr(kinship.matrix,"Var")
		kinship.matrix[upper.tri(kinship.matrix)] <- tmp[upper.tri(tmp)]
		ev <- eigen(kinship.matrix,symmetric=TRUE)$vectors
		rownames(ev) <- rownames(kinship.matrix)
		ev <- ev[data@phdata$id,1:naxes]
#		ev <- ev[tmeas,]
		rm(kinship.matrix,tmp);gc()
	} else {
		if (tmp[2]<naxes) stop("Eigenvector matrix supplied, but naxes supplied is smaller than naxes parameter")
		ev <- kinship.matrix[data@phdata$id,1:naxes]
#		ev <- ev[tmeas,]
		rm(kinship.matrix);gc()
	}
	resid <- lm(resid~ev)$resid

	lenn <- data@gtdata@nsnps;
	#out <- list()
	for (j in c(1:(times+1*(times>1)))) {
		if (j>1) resid <- sample(resid,replace=FALSE)
		chi2 <- .C("egscore",as.raw(data@gtdata@gtps),as.double(resid),as.integer(naxes),as.double(ev),as.integer(data@gtdata@nids),as.integer(data@gtdata@nsnps), as.integer(nstra), as.integer(strata), chi2 = double(7*data@gtdata@nsnps), PACKAGE="GenABEL")$chi2
#		if (any(data@gtdata@chromosome=="X")) {
#		  ogX <- data@gtdata[,data@gtdata@chromosome=="X"]
#		  sxstra <- strata; sxstra[ogX@male==1] <- strata[ogX@male==1]+nstra
#		  chi2X <- .C("egscore",as.raw(ogX@gtps),as.double(resid),as.integer(naxes),as.double(ev),as.integer(ogX@nids),as.integer(ogX@nsnps), as.integer(nstra*2), as.integer(sxstra), chi2 = double(7*ogX@nsnps), PACKAGE="GenABEL")$chi2
#		  revec <- (data@gtdata@chromosome=="X")
#		  revec <- rep(revec,6)
#		  chi2 <- replace(chi2,revec,chi2X)
#		  rm(ogX,chi2X,revec);gc(verbose=FALSE)
#		}
		if (j == 1) {
			chi2.1df <- chi2[1:lenn];
			###out$chi2.1df <- chi2.1df
			chi2.2df <- rep(NA,lenn) #chi2[(lenn+1):(2*lenn)];
			###out$chi2.2df <- chi2.2df
			actdf <- chi2[(2*lenn+1):(3*lenn)];
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
			effB <- chi2[(3*lenn+1):(lenn*4)]
			effAB <- chi2[(4*lenn+1):(lenn*5)]
			effBB <- chi2[(5*lenn+1):(lenn*6)]
			if (times>1) {
				pr.1df <- rep(0,lenn)
				pr.2df <- rep(0,lenn)
				pr.c1df <- rep(0,lenn)
			}
		} else {
			th1 <- max(chi2[1:lenn])
			pr.1df <- pr.1df + 1*(chi2.1df < th1)
			pr.2df <- pr.2df + 1*(chi2.2df < max(chi2[(lenn+1):(2*lenn)]))
			pr.c1df <- pr.c1df + 1*(chi2.c1df < th1)
			if (!quiet && ((j-1)/bcast == round((j-1)/bcast))) {
				cat("\b\b\b\b\b\b",round((100*(j-1)/times),digits=2),"%",sep="")
				flush.console()
			}
		}
	}
	if (times > bcast) cat("\n")

	if (times>1) {
		P1df <- pr.1df/times
		P1df <- replace(P1df,(P1df==0),1/(1+times))
		P2df <- pr.2df/times
		P2df <- replace(P2df,(P2df==0),1/(1+times))
		Pc1df <- pr.c1df/times
		Pc1df <- replace(Pc1df,(Pc1df==0),1/(1+times))
	} else {
		P1df <- pchisq(chi2.1df,1,lower.tail=FALSE)
		P2df <- pchisq(chi2.2df,actdf,lower.tail=FALSE)
		Pc1df <- pchisq(chi2.c1df,1,lower.tail=FALSE)
	}
	#out$lambda <- lambda
	#out$effB <- effB
	#out$effAB <- effAB
	#out$effBB <- effBB
	#out$snpnames <- data@gtdata@snpnames
	#out$map <- data@gtdata@map
	#out$chromosome <- data@gtdata@chromosome
	#out$idnames <- data@gtdata@idnames
	#out$formula <- match.call()
	#out$family <- paste("score test for association, eigen-adjustment")
	#out$N <- chi2[(6*lenn+1):(lenn*7)]
	#class(out) <- "scan.gwaa"
	#out
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
			family = "EIGENSCORETEST"
	) 
	out
}
