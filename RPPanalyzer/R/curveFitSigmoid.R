`curveFitSigmoid` <-
function(x,y, xVal=NULL, plot=F, detectionLimit=F)  {
	
    uniqueX <- sort(unique(x))

    if (detectionLimit)  {

        confLower = aggregate(y, by=list(x), function(xx) quantile(xx,0.025))[,2]
        confUpper = aggregate(y, by=list(x), function(xx) quantile(xx,0.975))[,2]
        med = aggregate(y, by=list(x), function(xx) median(xx))[,2]

        if(any(diff(med, 1) < 0))
            minX = max(min(uniqueX[confLower > confUpper[1]]), uniqueX[max(which(diff(med, 1) < 0))+1])
        else
            minX = min(uniqueX[confLower > confUpper[1]])

        if(minX == Inf)
            minX <- max(x)

        y <- y[x >= minX]
        x <- x[x >= minX]
    }

    if (length(unique(x)) > 1)  {

        #linear fit to find goos starting values
        helpFit <- try(rq(y ~ g.link(0.1*x)))

	if(inherits(helpFit, "try-error")) {
	    browser()
	}

        print(helpFit)

        #real fit
        #fit <- nlrq(y ~ curvePredictSigmoid(x, c(alpha, beta, gamma)),start=list(alpha=coef(helpFit)[1], beta=coef(helpFit)[2], gamma=0.1))
        #fit <- nlrq(y ~ curvePredictSigmoid(x, c(alpha, beta, gamma)))
        fit <- try(nlrq(y~curvePredictSigmoid(x, c(alpha,beta,gamma)), control=nlrq.control(maxiter=100, k=2, InitialStepSize=0.9, big=1e+20, eps=1e-07, beta=0.97), start=list(alpha=coef(helpFit)[1],beta=coef(helpFit)[2],gamma=0.1), trace=FALSE))

        #possible error ?
        if (inherits(fit,"try-error"))  {
            fit <- c(coef(helpFit, 0.1))
        }
        else {
            fit <- coef(fit)
        }

        if(is.null(xVal)) {
            val <- curvePredictSigmoid(max(x), fit)
        }
        else {
            val <- curvePredictSigmoid(xVal, fit)
        }

    }
    else {
        #no fit possible at all
        #returning the median of the the y-valuies, corresponding to x

        fit <- NULL
	#bic <- NULL

        if(is.null(xVal))
            xVal <- max(x)

        val <- median(y[x == xVal])

    }


    if (plot) {

        if (detectionLimit)  {
            #plot the lines indicating the detection limit
            abline(v=minX-0.1, col="red")
        }

        if(!is.null(fit)) {

            xPlot <- seq(min(x), max(x), 0.1)
            yPlot <- curvePredictSigmoid(xPlot, coef(fit))

            lines(xPlot, yPlot, col="green")

            #draw maximum of slope
            abline(v=-1/coef(fit)[3], col="blue")

        }



    }

    return(list(fit=fit, val=val))

}

