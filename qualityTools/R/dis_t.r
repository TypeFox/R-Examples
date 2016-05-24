setClass("distr", representation(x = "vector", name = "character", parameters = "numeric", sd = "numeric", n = "numeric", loglik = "numeric"))
setClass("distrCollection", representation(distr = "list"))
setMethod("[", signature(x = "distrCollection", i = "ANY"), function(x, i, drop = missing) {
    x@distr[i]
})
setMethod("show", signature(object = "distrCollection"), function(object) {
    distrList = object@distr
    cat("\n")
    for (i in seq(along = distrList)) {
        temp = distrList[[i]]
        cat("\n")
        cat("fitted distribution is", temp@name, ":\n")
        print(temp@parameters)
        cat("\n")
    }
})
setMethod("summary", signature(object = "distrCollection"), function(object) {
    numDist = length(object@distr)
    gofMatrix = data.frame(matrix(nrow = numDist, ncol = 3))
    names(gofMatrix) = c("Distribution", "A", "p.value")
    cat("\n------ Fitted Distribution and estimated parameters ------\n")
    for (i in seq(along = object@distr)) {
        distrObj = object@distr[[i]]
        x = distrObj@x
        distribution = distrObj@name
        parameters = distrObj@parameters
        statistic = NA
        p.value = NA
        temp = .myADTest(x, distribution)
        try(statistic <- as.numeric(temp$statistic), silent = TRUE)
        try(p.value <- as.numeric(temp$p.value), silent = TRUE)
        gofMatrix[i, ] = c(distribution, as.numeric(statistic), as.numeric(p.value))
        cat("\n")
        cat("fitted distribution is", distribution, ":\n")
        print(parameters)
    }
    cat("\n")
    cat("\n------ Goodness of Fit - Anderson Darling Test ------\n")
    cat("\n")
    gofMatrixPrint = gofMatrix
    gofMatrixPrint[, 2] = signif(as.numeric(gofMatrixPrint[, 2]), 4)
    gofMatrixPrint[, 3] = signif(as.numeric(gofMatrixPrint[, 3]), 4)
    print(gofMatrixPrint)
})
distribution = function(x, distribution = "weibull", start, ...) {
    #if (!require(MASS, quietly = TRUE)) 
     #   stop("Package MASS needs to be installed!")
    if (is.character(distribution)) 
        distribution = tolower(distribution)
    allDistr = c("beta", "cauchy", "chi-squared", "exponential", "f", "gamma", "geometric", "log-normal", "logistic", "negative binomial", "normal", "poisson", 
        "t", "weibull")
    if (distribution %in% allDistr) 
        distrVec = distribution
    else distrVec = c("normal")
    if (identical(distribution, "all")) 
        distrVec = allDistr
    if (identical(distribution, "quality")) 
        distrVec = c("normal", "log-normal", "exponential", "weibull")
    distrColl = new("distrCollection")
    for (i in seq(along = distrVec)) {
        fit = new("distr")
        temp = fitdistr(x, distrVec[i], ...)
        fit@x = x
        fit@name = distrVec[i]
        fit@parameters = temp$estimate
        fit@sd = temp$sd
        fit@loglik = temp$loglik
        distrColl@distr[i] = fit
    }
    return(distrColl)
}
setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
setMethod("plot", signature(x = "distr"), function(x, y, main, xlab, xlim, ylim, ylab, line.col, line.width, box = TRUE, ...) {
    object = x
    xVals = object@x
    parameters = object@parameters
    lq = NULL
    uq = NULL
    y = NULL
    if (missing(line.col)) 
        line.col = "red"
    if (missing(line.width)) 
        line.width = 1
    if (missing(main)) 
        main = object@name
    if (missing(xlab)) 
        xlab = "x"
    if (missing(ylab)) 
        ylab = "Density"
    distr = object@name
    qFun = .charToDistFunc(distr, type = "q")
    dFun = .charToDistFunc(distr, type = "d")
    adTestStats = .myADTest(xVals, distr)
    print(adTestStats)
    if (class(adTestStats) == "adtest") {
        A = adTestStats$statistic
        p = adTestStats$p.value
    }
    else {
        A = NA
        p = NA
    }
    histObj = hist(xVals, plot = FALSE)
    if (missing(xlim)) {
        lq = do.call(qFun, c(list(1e-04), as.list(parameters)))
        uq = do.call(qFun, c(list(0.9999), as.list(parameters)))
        xlim = range(lq, uq, xVals)
    }
    xPoints = seq(xlim[1], xlim[2], length = 200)
    yPoints = do.call(dFun, c(list(xPoints), as.list(parameters)))
    if (missing(ylim)) {
        ylim = range(0, histObj$density, yPoints)
    }
    hist(xVals, freq = FALSE, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, main = main, ...)
    lines(xPoints, yPoints, col = line.col, lwd = line.width)
    abline(h = 0)
    legend("topright", c(paste(c(names(parameters), "A", "p"), ": ", c(format(parameters, digits = 3), format(A, digits = 3), format(p, digits = 3))), sep = " "), 
        inset = 0.02)
    if (box) {
        box()
    }
})
.xylimits = function(distrCollection, lowerquantile = 0.001, upperquantile = 0.999) {
    x = NULL
    y = NULL
    for (i in seq(along = distrCollection@distr)) {
        object = distrCollection@distr[[i]]
        xValues = object@x
        parameters = object@parameters
        distr = object@name
        qFun = .charToDistFunc(distr, type = "q")
        dFun = .charToDistFunc(distr, type = "d")
        lq = do.call(qFun, c(list(lowerquantile), as.list(parameters)))
        uq = do.call(qFun, c(list(upperquantile), as.list(parameters)))
        x = range(xValues, x, lq, uq)
        histObj = hist(xValues, plot = FALSE)
        xPoints = seq(x[1], x[2], length = 200)
        yPoints = do.call(dFun, c(list(xPoints), as.list(parameters)))
        y = range(y, 0, histObj$density, yPoints)
    }
    invisible(list(xlim = x, ylim = y))
}
setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
setMethod("plot", signature(x = "distrCollection"), function(x, y, xlab, ylab, xlim, ylim, line.col, line.width, ...) {
    y = NULL
    object = x
    distrList = object@distr
    numDist = length(object@distr)
    numColWin = ceiling(numDist/2)
    if (missing(xlim)) 
        xlim = .xylimits(object)$xlim
    if (missing(ylim)) 
        ylim = .xylimits(object)$ylim
    if (missing(line.col)) 
        line.col = "red"
    if (missing(line.width)) 
        line.width = 1
    lapply(distrList, plot, xlim = xlim, ylim = ylim, line.col = line.col, line.width = line.width, ...)
    cat(paste("Total of", numDist, "plots created"))
    cat("\n")
    cat(paste("Use par(mfrow = c(2,", numColWin, ") to see all of them!", sep = ""))
    cat("\n")
})
qqPlot <- function (x, y, confbounds = TRUE, alpha, main, xlab, ylab, xlim, ylim, border = "red", bounds.col = "black", bounds.lty = 1,
    start, ...) 
{
    DB = FALSE
    parList = list(...)
    if (is.null(parList[["col"]])) 
        parList$col = 1:2
    if (is.null(parList[["pch"]])) 
        parList$pch = 19
    if (is.null(parList[["lwd"]])) 
        parList$lwd = 1
    if (is.null(parList[["cex"]])) 
        parList$cex = 1

    #if (!require(MASS)) 
    #    stop("Package MASS needs to be installed!")
    
    if (class(x) == "distrCollection") {
        distList = x@distr
        for (i in 1:length(distList)) {
            d = distList[[i]]
            do.call(qqPlot, c(list(x = d@x, y = d@name), parList))
        }
        invisible()
    }
    if (missing(y)) 
        y = "normal"
	if(missing(alpha))
		alpha = 0.05
    if (alpha <=0 || alpha >=1) 
            stop(paste("alpha should be between 0 and 1!"))		
    if (missing(main)) 
        main = paste("Q-Q Plot for", deparse(substitute(y)), 
            "distribution")
    if (missing(xlab)) 
        xlab = paste("Quantiles for", deparse(substitute(x)))
    if (missing(ylab)) 
        ylab = paste("Quantiles from", deparse(substitute(y)), 
            "distribution")
    if (is.numeric(y)) {
        cat("\ncalling (original) qqplot from namespace stats!\n")
        return(stats::qqplot(x, y, ...))
    }
    qFun = NULL
    theoretical.quantiles = NULL
    xs = sort(x)
    distribution = tolower(y)
    distWhichNeedParameters = c("weibull", "logistic", "gamma", 
        "exponential", "f", "geometric", "chi-squared", "negative binomial", 
        "poisson")
    
	
	# new
	threeParameterDistr = c("weibull3", "lognormal3", "gamma3")                
	threeParameter = distribution %in% threeParameterDistr
	if(threeParameter) distribution = substr(distribution, 1, nchar(distribution)-1)
	# end new

	if (is.character(distribution)) {
        qFun = .charToDistFunc(distribution, type = "q")
        if (is.null(qFun)) 
            stop(paste(deparse(substitute(y)), "distribution could not be found!"))
    }
    theoretical.probs = ppoints(xs)
	
    xq = NULL
    yq = quantile(xs, prob = c(0.25, 0.75))
    dots <- list(...)
    if (TRUE) {
        if (DB) 
            print("TODO: Pass the estimated parameters correctly")
        fitList = .lfkp(parList, formals(qFun))
        fitList$x = xs
        fitList$densfun = distribution
        if (!missing(start)) 
            fitList$start = start
        if (DB) {
            print(fitList)
            print("Ende")
        }
		# new
		if(!threeParameter){
        fittedDistr = do.call(fitdistr, fitList)
        parameter = fittedDistr$estimate

		#save the distribution parameter#
		thethas = fittedDistr$estimate
		# save the cariance-covariance matrix
		varmatrix = fittedDistr$vcov
		# end of my code
		
		# new code for three parameter
		} else {
			parameter = do.call(paste(".",distribution, "3", sep = ""), list(xs) )    ####
			threshold = parameter$threshold
		}
		
        parameter = .lfkp(as.list(parameter), formals(qFun))
        params = .lfkp(parList, formals(qFun))
        parameter = .lfrm(as.list(parameter), params)
        parameter = c(parameter, params)
        theoretical.quantiles = do.call(qFun, c(list(c(theoretical.probs)), 
            parameter))
		
		# new
		if(!threeParameter){		
		# array containing names of the distributions, for which conf intervals can be computed
		confIntCapable = c("exponential", "log-normal", "logistic", "normal", "weibull", "gamma", "beta", "cauchy")
		getConfIntFun = .charToDistFunc(distribution, type = ".confint")
		# if possible, compute the conf intervals
		if(confbounds == TRUE){
			if(distribution %in% confIntCapable){
				confInt = getConfIntFun(xs, thethas, varmatrix, alpha)
			}
		}# end of my code
		}
		
        xq <- do.call(qFun, c(list(c(0.25, 0.75)), parameter))
        if (DB) {
            print(paste("parameter: ", parameter))
            print(xq)
        }
    }
    else {
        params = .lfkp(parList, formals(qFun))
        params$p = theoretical.probs
        theoretical.quantiles = do.call(qFun, params)
        params$p = c(0.25, 0.75)
        xq = do.call(qFun, params)
    }

    params = .lfkp(parList, c(formals(plot.default), par()))	
	
	if(!threeParameter){
		params$y = theoretical.quantiles
	}  else {
		params$y = theoretical.quantiles+threshold
	}
	params$x = xs
    params$xlab = xlab
    params$ylab = ylab
    params$main = main
    if (!(is.null(params$col[1]) || is.na(params$col[1]))) 
        params$col = params$col[1]
    if (!missing(xlim)) 
        params$xlim = xlim
    if (!missing(ylim)) 
        params$ylim = ylim
    params$lwd = 1
    do.call(plot, params)
    pParams = params
    pParams = .lfkp(pParams, list(x = 1, y = 1, col = 1, cex = 1))
    do.call(points, pParams)
    params = .lfkp(parList, c(formals(abline), par()))
    params$a = 0
    params$b = 1
    params$col = border
    do.call(abline, params)
	
	if(!threeParameter){
	# plot the confInt if available
	if(confbounds == TRUE){
	 if(distribution %in% confIntCapable){
		params = .lfkp(parList, c(formals(lines), par()))	
		params$x = confInt[[3]]
		params$y = confInt[[1]]
		params$col = bounds.col
		params$lty = bounds.lty
		do.call(lines, params)
		
		params$x = confInt[[3]]
		params$y = confInt[[2]]
		params$col = bounds.col
		params$lty = bounds.lty
		do.call(lines, params)
		}
	} #end of my function
	}
	
    invisible(list(x = theoretical.quantiles, y = xs, int = params$a, 
        slope = params$b))
}

ppPlot <- function (x, distribution, confbounds = TRUE, alpha, probs, main, xlab, ylab, xlim, ylim, 
    border = "red", bounds.col = "black", bounds.lty = 1, grid = TRUE, box = TRUE, stats = TRUE, start, 
    ...) 
{
    DB = FALSE
    conf.level = 0.95
    conf.lines = TRUE
    #if (!require(MASS)) 
    #    stop("Package MASS needs to be installed!")
    if (!(is.numeric(x) | (class(x) == "distrCollection"))) 
        stop(paste(deparse(substitute(x)), " needs to be numeric or an object of class distrCollection"))
    parList = list(...)
    if (is.null(parList[["col"]])) 
        parList$col = c("black", "red", "gray")
    if (is.null(parList[["pch"]])) 
        parList$pch = 19
    if (is.null(parList[["lwd"]])) 
        parList$lwd = 1
    if (is.null(parList[["cex"]])) 
        parList$cex = 1
    if (DB) 
        print(parList)
    qFun = NULL
    xq = NULL
    yq = NULL
    x1 = NULL
	if(missing(alpha))
		alpha = 0.05
    if (alpha <=0 || alpha >=1) 
            stop(paste("alpha should be between 0 and 1!"))	
    if (missing(probs)) 
        probs = ppoints(11)
    else if (min(probs) <= 0 || max(probs) >= 1) 
        stop("probs should be values within (0,1)!")
    probs = round(probs, 2)
    if (is.numeric(x)) {
        x1 <- sort(na.omit(x))
        if (missing(xlim)) 
            xlim = c(min(x1) - 0.1 * diff(range(x1)), max(x1) + 
                0.1 * diff(range(x1)))
    }
    if (missing(distribution)) 
        distribution = "normal"
    if (missing(ylim)) 
        ylim = NULL
    if (missing(main)) 
        main = paste("Probability Plot for", deparse(substitute(distribution)), 
            "distribution")
    if (missing(xlab)) 
        xlab = deparse(substitute(x))
    if (missing(ylab)) 
        ylab = "Probability"
    if (class(x) == "distrCollection") {
        distList = x@distr
        for (i in 1:length(distList)) {
            d = distList[[i]]
            do.call(ppPlot, c(list(x = d@x, distribution = d@name), 
                parList))
        }
        invisible()
    }
    distWhichNeedParameters = c("weibull", "gamma", "logistic", 
        "exponential", "f", "geometric", "chi-squared", "negative binomial", 
        "poisson")
	# new
	threeParameterDistr = c("weibull3", "lognormal3", "gamma3")
	threeParameter = distribution %in% threeParameterDistr
	if(threeParameter) distribution = substr(distribution, 1, nchar(distribution)-1)
	# end new
	
    if (is.character(distribution)) {
        qFun = .charToDistFunc(distribution, type = "q")
        pFun = .charToDistFunc(distribution, type = "p")
        dFun = .charToDistFunc(distribution, type = "d")
        if (is.null(qFun)) 
            stop(paste(deparse(substitute(y)), "distribution could not be found!"))
    }
    dots <- list(...)
    if (TRUE) {
        if (DB) 
            print("TODO: Pass the estimated parameters correctly")
        fitList = .lfkp(parList, formals(qFun))
        fitList$x = x1
        fitList$densfun = distribution
        if (!missing(start)) 
            fitList$start = start
		
		if(!threeParameter){
        fittedDistr = do.call(fitdistr, fitList)
        parameter = fittedDistr$estimate
		#save the distribution parameter#
		thethas = fittedDistr$estimate
		# save the cariance-covariance matrix
		varmatrix = fittedDistr$vcov		
		} else {
			parameter = do.call(paste(".",distribution, "3", sep = ""), list(x1) )    ####
			print(parameter[3])
			threshold = parameter$threshold
		}
        parameter = .lfkp(as.list(parameter), formals(qFun))
        params = .lfkp(parList, formals(qFun))
        parameter = .lfrm(as.list(parameter), params)
        print(parameter)
        parameter = c(parameter, params)
        if (DB) {
            print(qFun)
            print(as.list(parameter))
            print(list(probs))
        }
		
				
		# new
		if(!threeParameter){		
		# array containing names of the distributions, for which conf intervals can be computed
		confIntCapable = c("exponential", "log-normal", "logistic", "normal", "weibull", "gamma", "beta", "cauchy")
		getConfIntFun = .charToDistFunc(distribution, type = ".confint")
		# if possible, compute the conf intervals
		if(confbounds == TRUE){
			if(distribution %in% confIntCapable){
				confInt = getConfIntFun(x1, thethas, varmatrix, alpha)
			}
		}# end of my code
		}
		

		y = do.call(qFun, c(list(ppoints(x1)), as.list(parameter)))
        yc = do.call(qFun, c(list(ppoints(x1)), as.list(parameter)))
        cv = do.call(dFun, c(list(yc), as.list(parameter)))
        print(cv)
        axisAtY = do.call(qFun, c(list(probs), as.list(parameter)))
        yq = do.call(qFun, c(list(c(0.25, 0.75)), as.list(parameter)))
        xq = quantile(x1, probs = c(0.25, 0.75))
        if (DB) {
            print(paste("parameter: ", parameter))
            print(xq)
        }
    }
    else {
        params = .lfkp(parList, formals(qFun))
        params$p = ppoints(x1)
        y = do.call(qFun, params)
        params$p = probs
        axisAtY = do.call(qFun, params)
        params$p = c(0.25, 0.75)
        yq = do.call(qFun, params)
        xq = quantile(x1, probs = c(0.25, 0.75))
    }
    params = .lfkp(parList, c(formals(plot.default), par()))

	params$x = x1
	params$y = y
    params$xlab = xlab
    params$ylab = ylab
    params$main = main
    params$xlim = xlim
    params$axes = FALSE
    params$lwd = 1
    if (!(is.null(params$col[1]) || is.na(params$col[1]))) 
        params$col = params$col[1]
    do.call(plot, params)
    pParams = params
    params = .lfkp(parList, list(cex.main = 1, cex.axis = 1, 
        cex.lab = 1))
    params$side = 1
    axisAtX = do.call(axis, params)
    params$side = 2
    params$at = axisAtY
    params$labels = probs
    params$las = 2
    do.call(axis, params)
    if (grid) {
        params = .lfkp(parList, c(formals(abline), list(lwd = 1, 
            col = 1)))
        params$h = axisAtY
        params$v = axisAtX
        if (!(is.null(params$col[3]) || is.na(params$col[3]))) 
            params$col = params$col[3]
        else params$col = 1
        if (!(is.null(params$lwd[2]) || is.na(params$lwd[2]))) 
            params$lwd = params$lwd[2]
        else params$lwd = 1
        do.call(abline, params)
    }
    pParams = .lfkp(pParams, list(x = 1, y = 1, col = 1, cex = 1))
    do.call(points, pParams)
    params = .lfkp(parList, c(formals(abline), par()))
    if(!threeParameter){
		params$a = 0
	}  else {
		params$a = -threshold
	}
    params$b = 1
    params$col = border
    do.call(abline, params)
	
	if(!threeParameter){
	# plot the confInt if available
	if(confbounds == TRUE){
	 if(distribution %in% confIntCapable){
		params = .lfkp(parList, c(formals(lines), par()))	
		params$x = confInt[[3]]
		params$y = confInt[[2]]
		params$col = bounds.col
		params$lty = bounds.lty
		do.call(lines, params)
		
		params$x = confInt[[3]]
		params$y = confInt[[1]]
		params$col = bounds.col
		params$lty = bounds.lty
		do.call(lines, params)		
		}
	} #end of my function
	}
    if (box) 
        box()
    invisible(list(x = x, y = y, int = params$a, slope = params$b))
}
