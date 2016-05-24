anoint.fit <- function(object,level=.05,interval=c(.5,3)){
	
	OBO <- obo(object)
	UIM <- uim(object)
	PIM.EXACT <- pim(object,exact=TRUE,interval=interval)	
	PIM.APPROX <- pim(object,exact=FALSE)
	
	#NUMBER ONE-BY-ONE
	K <- length(OBO$fit)
	
	#GLOBAL CRITERIA
	X2.95 <- qchisq((1-level),df=1)
	X2.975 <- qchisq((1-level/2),df=1)
	X2.adj <- qchisq((1-level/K),df=1)
	X2.975.adj <- qchisq((1-level/(2*K)),df=1)
	X2.95K <- qchisq((1-level),df=K)
	X2.975K <- qchisq((1-level/2),df=K)
	
	#REJECTIONS
	#ONE-STAGE
	obo.reject <- max(OBO$LRT)>=X2.95
	obo.adjust <- max(OBO$LRT)>=X2.adj
	uim.reject <- UIM$LRT>=X2.95K
	pim.exact.reject <- PIM.EXACT@LRT>=X2.95
	pim.approx.reject <- PIM.APPROX@LRT>=X2.95
	
	#TWO-STAGE
	pim.obo <- ifelse(PIM.EXACT@LRT>=X2.975,TRUE,
									max(OBO$LRT)>=X2.975.adj)
	
	pim.uim <- ifelse(PIM.EXACT@LRT>=X2.975,TRUE,
									UIM$LRT>=X2.975K)
	
	#P-VALUES
	obo.p <- pchisq(max(OBO$LRT),df=1,lower.tail=FALSE)
	uim.p <- pchisq(UIM$LRT,df=K,lower.tail=FALSE)
	pim.exact.p <- pchisq(PIM.EXACT@LRT,df=1,lower.tail=FALSE)
	pim.approx.p <- pchisq(PIM.APPROX@LRT,df=1,lower.tail=FALSE)
	
	#RESPONSE ESTIMATE
	theta.exact <- PIM.EXACT@coef$theta
	theta.approx <- PIM.APPROX@coef$theta
	
	new("anoint.fit",K=K,
		responsiveness = list(
			theta.exact = theta.exact,
			theta.approx = theta.approx),
		tests=list(
		obo.reject = obo.reject,
		obo.adjust = obo.adjust,
		uim.reject = uim.reject,
		pim.exact.reject = pim.exact.reject,
		pim.approx.reject = pim.approx.reject,
		pim.obo = pim.obo,
		pim.uim = pim.uim),
		pvalues=list(
		obo.p = obo.p,
		uim.p = uim.p,
		pim.exact.p = pim.exact.p,
		pim.approx.p = pim.approx.p
		),
		fits=list(
			obo = OBO,
			uim = UIM,
			pim.exact = PIM.EXACT,
			pim.approx = PIM.APPROX
		)
	)
}
