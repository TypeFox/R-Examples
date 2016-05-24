"mlreg.p" <- 
function(formula,data,snpsubset,idsubset,gtmode="additive",trait.type="guess")
#function(formula,data,snpsubset,idsubset,gtmode="additive",trait.type="guess",propPs=1.0)
{
	if (missing(data)) 
		stop("data argument can not be missing in lm.gwaa")
	if (!is(data,"gwaa.data")) 
		stop("data argument should have gwaa.data class")
	checkphengen(data)
	if (!missing(snpsubset)) data <- data[,snpsubset]
	if (!missing(idsubset)) data <- data[idsubset,]
	gtdata <- data@gtdata
	data <- data@phdata
#	if (!missing(snpsubset)) gtdata <- gtdata[,snpsubset]
#	if (!missing(idsubset)) {
#		gtdata <- gtdata[idsubset,]
#		data <- data[idsubset,]
#	}
	gc()
	
# from lm
	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
#	m <- match(c("formula", "data", "snpsubset", "idsubset"), names(mf), 0L)
	m <- match(c("formula", "data"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")
	x <- model.matrix(mt, mf)

# trait type (1=gaussian, 2=binomial, 3=survival)
	posttypes <- c("gaussian","binomial","survival")
	y <- model.extract(mf,"response")
	if (trait.type=="guess") {
		if (is(y,"numeric") || is(y,"integer")) {
			if (isbinomial(y)) ttype <- 2
			else ttype <- 1
		} else if (is(y,"GASurv")) {
			ttype <- 3
		} else {
			stop("cannot guess trait type: not numeric, not GASurv!")
		}
	} else {
		if (!any(posttypes==trait.type)) stop("trait type should be either \"gaussian\" or \"binomial\" or \"survival\"")
		ttype <- match(trait.type,posttypes)
	}
	if (ttype==2) test.type(y,"binomial")
	if (ttype==1) test.type(y,"gaussian")
# from lm.fit
	if (is.null(n <- nrow(x))) 
		stop("'x' must be a matrix")
	if (n == 0L) 
		stop("0 (non-NA) cases")
	p <- ncol(x)
	if (p == 0L)
		stop("no predictors in the model")
	ny <- NCOL(y)
	if (is.matrix(y) && ny == 1) 
		y <- drop(y)
	if (NROW(y) != n) 
		stop("incompatible dimensions")
	storage.mode(x) <- "double"
	storage.mode(y) <- "double"
# my code
# check mode 1=additive, 2=dominant 3=recessive 4=overdominant
	posgtmtypes <- c("additive","dominant","recessive","overdominant")
	if (!any(posgtmtypes==gtmode)) stop("gtmode should be one of: additive, dominant, recessive, overdominant")
	gtm <- match(gtmode,posgtmtypes)
# continue...
	gtdata <- gtdata[rownames(x),]
	if (gtdata@nids != n)
		stop("incompatible dimensions")
	if (any(gtdata@chromosome == "X") & !any(colnames(x)=="sex"))
		warning("You analysed X chromosome without adjustment for sex")
	if (ttype==1)
		chi2 <- .C("linreg_gwaa",as.double(y),as.double(x),as.raw(gtdata@gtps),as.integer(n),
			as.integer(ncol(x)),as.integer(gtdata@nsnps),as.integer(gtm),chi2=double(3*gtdata@nsnps))$chi2
	else if (ttype==2)
		chi2 <- .C("logreg_gwaa",as.double(y),as.double(x),as.raw(gtdata@gtps),as.integer(n),
			as.integer(ncol(x)),as.integer(gtdata@nsnps),as.integer(gtm),chi2=double(3*gtdata@nsnps))$chi2
	else if (ttype==3)
		chi2 <- .C("coxph_gwaa",as.double(y),as.double(x),as.raw(gtdata@gtps),as.integer(n),
			as.integer(ncol(x)),as.integer(gtdata@nsnps),as.integer(gtm),chi2=double(3*gtdata@nsnps))$chi2
	else 
		stop(paste("oppa, da oppa:",class(y),"!"))
	chi2 <- matrix(chi2,ncol=3)
	chi2 <- data.frame(chi2)
	colnames(chi2) <- c("N","beta","sebeta")
	rownames(chi2) <- gtdata@snpnames
	chi2[abs(chi2+999.9)<1.e-6] <- NA
	chi2.1df <- (chi2$beta/chi2$sebeta)^2
	out <- data.frame(snpnames=gtdata@snpnames,chromosome=gtdata@chromosome,map=gtdata@map,
			N=chi2$N,effB=chi2$beta,se.effB=chi2$sebeta,chi2.1df=chi2.1df,
			P1df=pchisq(chi2.1df,1,lower.tail=F),stringsAsFactors = FALSE)
#	class(out) <- "scan.gwaa"
	out
}

