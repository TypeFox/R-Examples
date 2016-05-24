`curveFitLinear` <-
function(x, y, xVal=NULL, method="quantreg", plot=F, detectionLimit=F)  {

    uniqueX <- sort(unique(x))

    if (detectionLimit)  {

	#use the quantiles 0.025 and 0.975 to determine a lower and an upper bound
	#for each dilutions steps
	# 1. criterion: lower bound of the first measurable step has to be higer than the upper bound of the last
	#               step (is considered pur background noise)
	# 2. criterion: medians of measurable dilution steps have to be rising
        confLower = aggregate(y, by=list(x), function(xx) quantile(xx,0.025))[,2]
        confUpper = aggregate(y, by=list(x), function(xx) quantile(xx,0.975))[,2]
        med = aggregate(y, by=list(x), function(xx) median(xx))[,2]

        if(any(diff(med, 1) < 0))  {
            minX = max(min(uniqueX[confLower > confUpper[1]]), uniqueX[max(which(diff(med, 1) < 0))+1])
        }
        else  {
            minX = min(uniqueX[confLower > confUpper[1]])
        }

        if(minX == Inf)  {
            minX <- max(x)
        }

        y <- y[x >= minX]
        x <- x[x >= minX]
    }


    if (method=="quantreg") {

        fit <- try(rq(y~x))
        if ( inherits(fit,"try-error")){
       fit=NULL
        }
    }
    else if (method=="simple")  {

        fit <- try(lm(y~x))
        if ( inherits(fit,"try-error")){
       fit=NULL
        }
    }

    if (!is.null(fit)){

    bic <- AIC(fit, k=log(length(x)))

    print(fit)

    fit <- coef(fit)

    if(is.null(xVal))
        val <- curvePredictLinear(max(x), fit)
    else
        val <- curvePredictLinear(xVal, fit)
    }
    else{
    bic=NULL
    val=median(y[x==max(x)])
    }

    if(plot)  {

        if(detectionLimit)
            abline(v=minX-0.1, col="red")

        if(!is.null(fit))
            lines(x=range(x), y=curvePredictLinear(x=range(x), params=coef(fit)) ,col="green")

    }


    return(list(fit=fit, val=val, bic=bic))


}

