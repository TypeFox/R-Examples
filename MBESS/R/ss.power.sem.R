ss.power.sem <- function (F.ML=NULL, df=NULL, RMSEA.null=NULL, RMSEA.true=NULL, F.full=NULL, F.res=NULL, RMSEA.full=NULL, RMSEA.res=NULL, 
	df.full=NULL, df.res=NULL, alpha=.05, power=.80) {
	
    if (power < 0 | power > 1) 
        stop("The value of 'power' must be between 0 and 1")

	if (alpha>1 | alpha<0) stop ("The value of 'alpha' must be between 0 and 1")
	
	input.F.ml <- input.RMSEA.true <- input.F.full <- input.RMSEA.full <- 0
	One.model <- Two.model <- FALSE
	
	if (!is.null(F.ML)) {
		if(is.null(df)) stop ("Because 'F.ML' is specified, the function requires the model's degrees of freedom 'df' as input as well.")
		rmsea0 <- 0
		rmseaa <- sqrt(F.ML/df)
		d <- df
		input.F.ml <- 1
		One.model <- TRUE
		} else input.F.ml <- 0
	
	if (!is.null(RMSEA.true)) {
		if(is.null(RMSEA.null) || is.null(df)) stop ("Because 'RMSEA.true' is specified, the function requires 'RMSEA.null' and 'df' as input as well.")
		rmsea0 <- RMSEA.null
		rmseaa <- RMSEA.true
		d <- df
		input.RMSEA.true <-1
		One.model <- TRUE
		} else input.RMSEA.true <- 0
	
	if(!is.null(F.full)){
		if(is.null(df.full) || is.null(df.res) || is.null(F.res) ) stop ("Because 'F.full' is specified, the function requires 'F.res', 'df.full', and 'df.res' as input as well")
		if(!is.null(F.ML) || !is.null(RMSEA.full) || !is.null(RMSEA.true) ) stop("Because 'F.full' is specified, do not put in 'F.ML' or any RMSEA values")
		rmseaa <- sqrt(F.res/df.res)
		rmseab <- sqrt(F.full/df.full)
		da <- df.res
		db <- df.full
		input.F.full <-1
		Two.model <- TRUE
		} else input.F.full <- 0 
		
	if(!is.null(RMSEA.full)){
		if(is.null(RMSEA.res) || is.null(df.full) || is.null(df.res)) stop("Because 'RMSEA.full' is specified, the function also requires 'RMSEA.res', 'df.full', and 'df.res' as input")
		if(RMSEA.full > RMSEA.res) stop ("The full model's RMSEA must be smaller than or equal to the restricted model's RMSEA")
		rmseaa <- RMSEA.res
		rmseab <- RMSEA.full
		da <- df.res
		db <- df.full
		input.RMSEA.full <-1
		Two.model <- TRUE
		} else input.RMSEA.full <- 0	
		
	if (input.F.ml+input.RMSEA.true+input.F.full+input.RMSEA.full==0) stop ("Please specify at least")	
	
	alpha <- alpha
	desired <- power
	
	if(One.model){	
		pow <- 0.0
		n <- 0

		while (pow<desired) {
			n <- n+100
			ncp0 <- (n-1)*d*rmsea0^2
			ncpa <- (n-1)*d*rmseaa^2
  
		if(rmsea0<rmseaa) {
			cval <- qchisq(alpha,d,ncp=ncp0,lower.tail=F)
			pow <- pchisq(cval,d,ncp=ncpa,lower.tail=F)
			}
		else {
			cval <- qchisq(1-alpha,d,ncp=ncp0,lower.tail=F)
			pow <- 1-pchisq(cval,d,ncp=ncpa,lower.tail=F)
			}
		}


		foo <- -1
		newn <- n
		interval <- 200
		powdiff <- pow - desired

		while (powdiff>.001) {
			interval <- interval*.5
			newn <- newn + foo*interval*.5
			ncp0 <- (newn-1)*d*rmsea0^2
			ncpa <- (newn-1)*d*rmseaa^2
  
			if(rmsea0<rmseaa) {
				cval <- qchisq(alpha,d,ncp=ncp0,lower.tail=F)
				pow <- pchisq(cval,d,ncp=ncpa,lower.tail=F)
				}
			else {
				cval <- qchisq(1-alpha,d,ncp=ncp0,lower.tail=F)
				pow <- 1-pchisq(cval,d,ncp=ncpa,lower.tail=F)
				}
			powdiff <- abs(pow-desired)
			if (pow<desired) {foo <- 1}
			if (pow>desired) {foo <- -1}
		}

	minn <- newn
	}
	
###########################################
	if(Two.model){
		ddiff <- da-db
		fa <- da*rmseaa^2
		fb <- db*rmseab^2

		pow <- 0.0
		n <- 0

		while (pow<desired) {
			n <- n+100
			ncp <- (n-1)*(fa-fb)
		   
			cval <- qchisq(alpha,ddiff,ncp=0,lower.tail=F)
			pow <- pchisq(cval,ddiff,ncp=ncp, lower.tail=F)
		}


		foo <- -1
		newn <- n
		interval <- 200
		powdiff <- pow - desired
		while (powdiff>.001) {
			interval <- interval*.5
			newn <- newn + foo*interval*.5
			ncp <- (newn-1)*(fa-fb)
			#compute power
			cval <- qchisq(alpha,ddiff,ncp=0,lower.tail=F)
			pow <- pchisq(cval,ddiff,ncp=ncp,lower.tail=F)

			if(pow<desired) {foo <- 1}
			else {foo <- -1}

			powdiff <- abs(pow-desired)
		}
		minn <- newn
	}
	
return(ceiling(minn))
}