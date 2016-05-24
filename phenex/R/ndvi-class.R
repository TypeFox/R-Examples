.initNDVI <-function() {
	setClass("OptionalInteger")
	setClassUnion("OptionalInteger", c("logical", "integer"))

	setClass("OptionalVector")
	setClassUnion("OptionalVector", c("NULL", "vector"))

	setClass("NDVI", representation(year="OptionalInteger", values="vector", correctedValues="OptionalVector", modelledValues="OptionalVector"), 
		prototype(year=NA, values=rep(NA, 365), correctedValues=NULL, modelledValues=NULL))

	setGeneric("year", function(x) standardGeneric("year"))
	setMethod("year", "NDVI", function(x) x@year)	

	setGeneric("values", function(x) standardGeneric("values"))
	setMethod("values", "NDVI", function(x) x@values)

	setGeneric("correctedValues", function(x) standardGeneric("correctedValues"))
	setMethod("correctedValues", "NDVI", function(x) x@correctedValues)

	setGeneric("modelledValues", function(x) standardGeneric("modelledValues"))
	setMethod("modelledValues", "NDVI", function(x) x@modelledValues)

	setGeneric("isLeapYear", function(x) standardGeneric("isLeapYear"))
	setMethod("isLeapYear", "NDVI", function(x) {
		leap <- FALSE
		y <- year(x)
		if ((y %% 4)==0) { leap <- TRUE }
		if ((y %% 100)==0) { leap <- FALSE }
		if ((y %% 400)==0) { leap <- TRUE }
		return(leap)
	})

	setValidity("NDVI", function(object) {
		if (!is.na(year(object)) & !is.integer(year(object))){
			stop("Year has to be 'integer' or 'NA'")
		}
		if (!is.vector(values(object))){
			stop("'values' has to be a vector")
		}
		if (length(which(!is.na(values(object)))) > 0){
			if ( (min(values(object), na.rm=TRUE) < -1) | 
				(max(values(object), na.rm=TRUE) > 1)){
				stop("Values of 'values' have to be in interval (-1,1)")
			}
		}
		if ( !is.null(correctedValues(object)) & !is.vector(correctedValues(object)) ){
			stop("'correctedValues' has to be NULL or a vector")
		}
		if ( !is.null(modelledValues(object)) & !is.vector(modelledValues(object)) ){
			stop("'correctedValues' has to be NULL or a vector")
		}
		return(TRUE)
	})

	setGeneric("year<-", function(x, value) standardGeneric("year<-"))
	setReplaceMethod("year", "NDVI", function(x, value) {x@year <- value; validObject(x); x})

	setGeneric("values<-", function(x, value) standardGeneric("values<-"))
	setReplaceMethod("values", "NDVI", function(x, value) {x@values <- value; validObject(x); x})

	setGeneric("runningAvg<-", function(x, value) standardGeneric("runningAvg<-"))
	setReplaceMethod("runningAvg", "NDVI", function(x, value) {x@runningAvg <- value; validObject(x); x})

	setGeneric("correctedValues<-", function(x, value) standardGeneric("correctedValues<-"))
	setReplaceMethod("correctedValues", "NDVI", function(x, value) {x@correctedValues <- value; validObject(x); x})

	setGeneric("modelledValues<-", function(x, value) standardGeneric("modelledValues<-"))
	setReplaceMethod("modelledValues", "NDVI", function(x, value) {x@modelledValues <- value; validObject(x); x})

	setGeneric("bise1", function(x, slidingperiod) standardGeneric("bise1"))
	setMethod("bise1", "NDVI", function(x, slidingperiod) {
		if (is.na(year(x))){ stop("'year' has to be set") }
		if (missing("slidingperiod")){ slidingperiod <- 40 }
		days <- ifelse(isLeapYear(x), 366, 365)

		res <- correctedndvi <- vector(mode="numeric",length=days);

		ndvi <- values(x)
		ndvi <- ifelse( is.na(ndvi), 0, ndvi)

		res <- .C("bise", rdays=as.integer(days), ndvi=as.numeric(ndvi), 
			rslidperiod=as.integer(slidingperiod), cndvi=as.numeric(correctedndvi), 
			PACKAGE="phenex")$cndvi

		res <- ifelse ( res <= 0, NA, res)

		#evaluate res (not enough values)
		if (length(which(is.na(res) == FALSE)) < 5){
			correctedValues(x) <- as.vector(rep(NA, days))
			validObject(x)
			return(x)
		}
		if (res[order(res, decreasing=TRUE)[5]] < 0.1){
			correctedValues(x) <- as.vector(rep(NA, days))
			validObject(x)
			return(x)
		}

		peaks <- c()
		threshold <- 1.5
		check <- which(is.na(res) == FALSE)
		for (i in 1:length(check)){
			if (i==1) {
				if (res[check[i]] > (threshold * mean(res[c(check[length(check)], check[i+1])]))){
					if (res[check[i]] > 0.3) {
						peaks <- c(peaks,check[i])
					}
				}
			} else {
				if (i==length(check)){
					if (res[check[i]] > (threshold * mean(res[c(check[i-1], check[1])]))){
						if (res[check[i]] > 0.3) {
							peaks <- c(peaks,check[i])
						}
					}
				} else {
					if (res[check[i]] > (threshold * mean(res[c(check[i-1], check[i+1])]))){
						if (res[check[i]] > 0.3) {
							peaks <- c(peaks,check[i])
						}
					}
				}
			}
		}
		res[peaks] <- NA
		correctedValues(x) <- as.vector(res)
		validObject(x)
		return(x)
	})

	setGeneric("bise2", function(x, slidingperiod, growthFactorThreshold, cycleValues) standardGeneric("bise2"))
	setMethod("bise2", "NDVI", function(x, slidingperiod, growthFactorThreshold, cycleValues) {
		if (is.na(year(x))){ stop("'year' has to be set") }
		if (missing("slidingperiod")){ slidingperiod <- 40 }
		if (missing("growthFactorThreshold")){ growthFactorThreshold <- 0.1 }
		if (missing("cycleValues")){ cycleValues <- TRUE }
		daysofyear <- ifelse(isLeapYear(x), 366, 365)

		ndvi <- values(x)
		ndvi[which(ndvi <= 0)] <- NA

		if (cycleValues){ 
			days <- 3*daysofyear
			ndvi <- c(ndvi,ndvi,ndvi)
		} else {
			days <- daysofyear
		}
		
		corrected <- rep(NA,length=days)
		pos <- 1
		lastValidPos <- 0
		while (pos <= days){
			if (is.na(ndvi[pos])){
				pos <- pos+1
				next
			}

			if (lastValidPos==0){
				if ((ndvi[pos] > mean(ndvi, na.rm=TRUE)/5) && 
					(ndvi[pos] <= mean(ndvi, na.rm=TRUE))){
					corrected[pos] <- ndvi[pos]
					lastValidPos <- pos
				}
				pos <- pos+1
				next
			}

			validIncrease <- (1+(growthFactorThreshold*(pos-lastValidPos)))*
						corrected[lastValidPos]

			validIncrease <- ifelse(validIncrease > 1, 1, validIncrease)

			if (ndvi[pos] >= corrected[lastValidPos]) {
				if ((ndvi[pos] <= validIncrease)) {
					corrected[pos] <- ndvi[pos]
					lastValidPos <- pos
				}
				pos <- pos+1
			} else if (ndvi[pos] < corrected[lastValidPos]) {

				endOfPeriod <- pos+slidingperiod
				endOfPeriod <- ifelse(endOfPeriod>days,days,endOfPeriod)
				period <- pos:endOfPeriod
				nextValues <- ndvi[period]

				cloudValues <- nextValues-corrected[lastValidPos]
				cloudHop <- which(cloudValues>0)
				if (length(cloudHop)>0){
					pos <- period[cloudHop[1]]
					next
				}

				slopeThreshold <- ndvi[pos]+0.2*(corrected[lastValidPos]-ndvi[pos])
				values <- nextValues-slopeThreshold
				possibleValues <- which(values>0)
				if (length(possibleValues)>0){
					pos <- period[possibleValues[1]]
				} else {
					corrected[pos] <- ndvi[pos]
					lastValidPos <- pos
					pos <- pos+1
				}
			}
		}

		if (cycleValues){
			corrected <- corrected[(daysofyear+1):(2*daysofyear)]
		}

		correctedValues(x) <- corrected
		validObject(x)
		return(x)
	})

	setGeneric("bise", function(x, slidingperiod, growthFactorThreshold) standardGeneric("bise"))
	setMethod("bise", "NDVI", function(x, slidingperiod, growthFactorThreshold) {
		if (is.na(year(x))){ stop("'year' has to be set") }
		if (missing("slidingperiod")){ slidingperiod <- 40 }
		if (missing("growthFactorThreshold")){ growthFactorThreshold <- 0.1 }

		x <- bise2(x, slidingperiod, growthFactorThreshold)

		validObject(x)
		return(x)
	})

	setGeneric("runningAvg", function(x, window) standardGeneric("runningAvg"))
	setMethod("runningAvg", "NDVI", function(x, window) {
		if (is.na(year(x))){ stop("'year' has to be set") }
		if (missing("window")){ window <- 7 }
		days <- ifelse(isLeapYear(x), 366, 365)

		res <- correctedndvi <- vector(mode="numeric",length=days);

		ndvi <- values(x)
		ndvi <- ifelse( is.na(ndvi), -1, ndvi)

		cndvi <- .C("runAVG", rdays = as.integer(days), 
				ndvi=as.numeric(ifelse(is.na(ndvi), -1, ndvi)),
				window=as.integer(7), cndvi=as.numeric(correctedndvi), 
				PACKAGE="phenex")$cndvi
		res <- ifelse(cndvi < 0, NA, cndvi)

		#evaluate res (not enough values)
		if (length(which(is.na(res) == FALSE)) < 5){
			correctedValues(x) <- as.vector(rep(NA, days))
			validObject(x)
			return(x)
		}
		if (res[order(res, decreasing=TRUE)[5]] < 0.1){
			correctedValues(x) <- as.vector(rep(NA, days))
			validObject(x)
			return(x)
		}

		correctedValues(x) <- as.vector(res)
		validObject(x)
		return(x)
	})

	setGeneric("modelValues", function(x, method, ...) standardGeneric("modelValues"))
	setMethod("modelValues", "NDVI", function(x, method, ...) {
		# ...: filter.threshold for fft
		# ...: asym for gauss
		# ...: window.sav, degree, smoothing for savgol
		# method can be 'LinIP', 'Spline', 'DSig', 'DSigC',
		# 'DLogistic', 'Gauss', 'Growth', 'FFT' or 'SavGol'
		if (tolower(method)!="linip" & tolower(method)!="spline" &
		tolower(method)!="dsig" & tolower(method)!="dsigc" &
		tolower(method)!="dlogistic" & tolower(method)!="gauss" &
		tolower(method)!="growth" & tolower(method)!="fft" &
		tolower(method)!="savgol") {
			stop("method has to be one of the following: 
				'LinIP', 'Spline', 'DSig', 'DSigC', 
				'DLogistic', 'Gauss', 'Growth', 
				'FFT' or 'SavGol'")
		}
		arguments <- names(list(...))
		if (tolower(method)=="gauss"){
			if (is.na(match("asym",arguments))){
				asym <- FALSE
			} else { asym <- list(...)$asym }
		}
		if (tolower(method)=="fft"){
			if (is.na(match("filter.threshold",arguments))){
				filter.threshold <- 3
			} else { filter.threshold <- list(...)$filter.threshold }
		}
		if (tolower(method)=="savgol"){
			if (is.na(match("window.sav",arguments))){
				window.sav <- 7
			} else { window.sav <- list(...)$window.sav }
			if (is.na(match("degree",arguments))){
				degree <- 2
			} else { degree <- list(...)$degree }
			if (is.na(match("smoothing",arguments))){
				smoothing <- 10
			} else { smoothing <- list(...)$smoothing }
		}

		if (is.na(year(x))){ stop("'year' has to be set") }
		
		ndvi.mod <- rep(NA, ifelse(isLeapYear(x),366,365))
		if(is.null(correctedValues(x))){
			ndvi <- values(x)
		} else {
			ndvi <- correctedValues(x)
		}
		
		if (length(which(!is.na(ndvi))) >= 5) {
			if (tolower(method)=="linip"){ ndvi.mod <- .linIP(ndvi) }
			if (tolower(method)=="spline"){ ndvi.mod <- .spline(ndvi) }
			if (tolower(method)=="dsig"){ ndvi.mod <- .dSig(ndvi) }
			if (tolower(method)=="dsigc"){ ndvi.mod <-  .dSigC(ndvi) }
			if (tolower(method)=="dlogistic"){ ndvi.mod <- .dLogistic(ndvi) }
			if (tolower(method)=="gauss"){ ndvi.mod <- .gauss(ndvi, asym) }
			if (tolower(method)=="growth"){ ndvi.mod <- .growth(ndvi) }
			if (tolower(method)=="fft"){ ndvi.mod <- .fftfilter(ndvi, filter.threshold) }
			if (tolower(method)=="savgol"){ ndvi.mod <- .savGol(ndvi, window.sav, degree, smoothing) }
		}
		modelledValues(x) <- ndvi.mod
		validObject(x)
		return(x)
	})

	if (!isGeneric("plot"))
      		setGeneric("plot", function(x,y,...) standardGeneric("plot"))
		setMethod("plot", "NDVI", function(x,y=NULL,...) {
		plot(1:ifelse(isLeapYear(x),366,365), values(x), type="p",
			xlim=c(0,365), ylim=c(0,1), xlab="Day of the Year",
			ylab="NDVI", col="black", ...)
		if(!is.null(correctedValues(x))){
			points(1:ifelse(isLeapYear(x),366,365), correctedValues(x), 
				col="red", ...)
		}
		if (!is.null(modelledValues(x))){
			lines(1:ifelse(isLeapYear(x),366,365), modelledValues(x), 
				col="blue", ...)
		}
	})

	setGeneric("phenoPhase", function(x, phase, method, threshold) standardGeneric("phenoPhase"))
	setMethod("phenoPhase", "NDVI", function(x, phase, method, threshold) {
		if ( length(which(!is.na(modelledValues(x)))) < 5 ) { return(NA) }
		if (tolower(phase)!="max" & tolower(phase)!="min" & 
			tolower(phase)!="greenup" & tolower(phase)!="senescence"){
			stop("'phase' has to be 'min', 'max', 'greenup' or 'senescence'")
		}
		
		if (tolower(phase)=="max"){
			return(order(modelledValues(x), decreasing=TRUE, na.last=TRUE)[1])
		}
		if (tolower(phase)=="min"){
			maxDOY <- order(modelledValues(x), decreasing=TRUE, na.last=TRUE)[1]
			return(order(modelledValues(x)[1:maxDOY], decreasing=FALSE, na.last=TRUE)[1])
		}

		if (tolower(method)!="local" & tolower(method)!="global") {
			stop("'method' has to be 'local' or 'global'")
		}
		if (!is.numeric(threshold)){
			stop("'threshold' has to be a numeric value")
		}
		if (threshold < 0 | threshold > 1){
			stop("'threshold' has to be in interval (0,1)")
		}

		if (tolower(phase)=="greenup"){ start <- 1 }
		if (tolower(phase)=="senescence"){ start <- length(modelledValues(x)) }

		doy <- NA
		model <- modelledValues(x) 

		if (method=="local"){
			maxDoy <- order(model, decreasing=TRUE, na.last=TRUE)[1]
			if (tolower(phase)=="senescence"){
				minmodel <- min(model[maxDoy:length(model)], na.rm=TRUE)
			} else {
				minmodel <- min(model[1:maxDoy], na.rm=TRUE)
			}
			maxmodel <- max(model, na.rm=TRUE)
			threshold <- ((maxmodel - minmodel) * threshold) + minmodel
		}

		modelhalf <- model[start:which(model == max(model, na.rm=TRUE))[1]]
		thresvec <- modelhalf-rep(threshold, length(modelhalf))

		if (length(which(thresvec > 0)) > 0){
			doy <- which(thresvec > 0)[na.omit(match(which(thresvec <= 0)+1, which(thresvec > 0)))]
			doy <- doy[order(doy, decreasing=FALSE)[1]]
	
			if (tolower(phase)=="senescence"){ doy <- length(modelledValues(x)) + 1 - doy}
		} else {
			doy <- NA
		}

		return(doy)
	})

	
	setGeneric("rsquare", function(x) standardGeneric("rsquare"))
	setMethod("rsquare", "NDVI", function(x) {
		if (is.null(modelledValues(x))){ return(NA) }
		if(is.null(correctedValues(x)) | length(correctedValues(x)) < 2){
			ndvi <- values(x)
		} else {
			ndvi <- correctedValues(x)
		}

		model <- modelledValues(x)

		check <- which(!is.na(ndvi) & !is.na(model))

		if ((length(check) > 0) && 
			(var(ndvi[check])!=0) && 
			(var(model[check])!=0)){
			r2 <- cor(ndvi[check], model[check])^2
		} else {
			r2 <- NA
		}
		return(r2)
	})

	setGeneric("integrateTimeserie", function(x, start, end) standardGeneric("integrateTimeserie"))
	setMethod("integrateTimeserie", "NDVI", function(x, start, end){
		errorreturn <- list(value=NA, abs.error=NA, subdivisions=NA, message=NA, call=NA)
		if(is.null(modelledValues(x))){ return(errorreturn) }
		if (is.na(start) | is.na(end)){ return(errorreturn) }
		model <- modelledValues(x)
		if(length(which(!is.na(model)))<2){ return(errorreturn) }
		modelfunc <- approxfun(x=(-length(model)+1):(2*length(model)), y=c(model,model,model))
		count <- 1
		repeat {
			subdiv <- 100
			res <- try(integrate(f=modelfunc, lower=start, 
					upper=end, subdivisions=subdiv), 
				silent=TRUE)
			if (!inherits(res, "try-error")){
				return(res)
			} else {
				subdiv <- subdiv*2
				count <- count+1
				if (count >= 20){
					break;
				}
			}
		}
		return(errorreturn)
	})
}
