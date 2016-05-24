
##' @title arma.forecast
##' 
##' @description Forecasting of (multivariate) time series of
##' marima type using marima type model.
##'
##' @param series = matrix holding the kvar-variate timeseries.
##' The series is assumed to have the same format
##' as the timeseries analysed by marima BEFORE differencing (if
##' differencing was used via define.dif)
##' (the length, though, does not need to be the same but can be shorter
##' or longer). Results
##' from estimating the model (for the differenced data, if used) are assumed
##' to be saved in the input-object 'marima' (see 'usage') by marima.
##'
##' The series is assumed to have the total length=(nstart+nstep) (but it
##' may be longer. In any case the forecasting is starting from nstart
##' continuing to nstart+nstep. Future values already present or initialised,
##' for example, as NAs are overwritten with the forecasted values.)
##' 
##' An example of a series prepared for forcasting is in the marima library:
##' 'data(austr)': (see below, the example).
##'
##' If future (independent) x-values for the forecasting are to be used
##' these values must be supplied in 'series' at the proper places before
##' calling 'arma.forecast(...)' (that is except the x-value(s)
##' corresponding to the last prediction). 
##'
##' @param marima = the object holding the marima results to be used for
##' the forecasting, that is an output object created by marima.
##'
##' If the ar- and/or the ma-model do not include a leading unity matrix
##' this is automatically taken care of in the function (in that case the
##' dimensions of the model arrays used will be, respectively,
##' (kvar,kvar,p+1) and (kvar,kvar,q+1)) after inserting the leading
##' unity matrix (if the object 'marima' was produced by marima, this
##' will automatically be OK.
##'
##' @param nstart = starting point for forecasting (1st forecast values
##' will be for time point t = nstart+1).
##'
##' @param nstep = length of forecast (forecasts will be for time points
##' nstart+1,...,nstart+nstep).
##'
##' @param dif.poly (most often) output from the function define.dif holding
##' the ar-representation of the differencing polynomial
##' (define.dif$dif.poly).
##' If a differenced timeseries was analysed by marima
##' the forecast-variance/covariance matrices are calculated for the
##' aggregated (original) timeseries if 'dif.poly' is specified. If not,
##' the forecast-variance/covariance matrices are calculated for the
##' differenced time series. If forecasting is wanted for the original
##' (not differenced) time series the 'dif.poly' created by define.dif
##' must be specified.
##'
##' @return forecasts = forecasted values following the nstart first values
##' of the input series (at time points 'nstart+1,...,nstart+nstep').
##' The forecasted values will be (over-) written in the input series at
##' the proper future positions (if relevant).
##'
##' @return residuals = corresponding residuals for input series followed by
##' nstep future residuals (all=0).
##'
##' @return prediction.variances = (kvar,kvar,nstep) array containing
##' prediction covariance matrices corresponding to the nstep forecasts.
##' @return nstart = starting point for prediction (1st prediction at point
##' nstart+1).
##' @return nstep = length of forecast
##'
##' @examples
##'
##' library(marima)
##' data(austr)
##' series<-t(austr)
##' Model5 <- define.model(kvar=7, ar=1, ma=1, rem.var=1, reg.var=6:7)
##' Marima5 <- marima(series[, 1:90], Model5$ar.pattern, Model5$ma.pattern, 
##' penalty=1)
##' nstart  <- 90
##' nstep   <- 10
##' Forecasts <- arma.forecast(series=series, marima=Marima5, 
##'                nstart=nstart, nstep=nstep)
##' Year<-series[1, 91:100];
##' Predict<-Forecasts$forecasts[2, 91:100]
##' stdv<-sqrt(Forecasts$pred.var[2, 2, ])
##' upper.lim=Predict+stdv*1.645
##' lower.lim=Predict-stdv*1.645
##' Out<-rbind(Year, Predict, upper.lim, lower.lim)
##' print(Out)
##' # plot results:
##' plot(series[1, 1:100], Forecasts$forecasts[2, ], type='l', xlab='Year', 
##' ylab='Rate of armed suicides', main='Prediction of suicides by firearms', 
##' ylim=c(0.0, 4.1))
##' lines(series[1, 1:90], series[2, 1:90], type='p')
##' grid(lty=2, lwd=1, col='black')
##' Years<-2005:2014
##' lines(Years, Predict, type='l')
##' lines(Years, upper.lim, type='l')
##' lines(Years, lower.lim, type='l')
##' lines(c(2004.5, 2004.5), c(0.0, 2.0), lty = 2)
##'
##' @importFrom graphics grid plot
##' 
##' @export

arma.forecast <- function( series=NULL, marima=NULL, 
                           nstart=NULL, nstep=1, dif.poly=NULL ) {
    Y <- series[, c( 1:(nstart + nstep) )]
    d <- dim(Y)
    if (d[1] > d[2]) {
        Y <- t(Y)
        d <- dim(Y)
    }

    if (is.null(dif.poly)) {
        dif.poly <- array(diag(d[1]), dim = c( d[1], d[1], 1 ))
    }

    dif.poly <- check.one( dif.poly )

    colnames(Y) <- c(1:d[2])

    means    <- marima$mean.pattern
    averages <- marima$averages
    Constant <- marima$Constant

    if (is.null(nstart)) {
        nstart <- d[2]
    }
    
    Names <- rownames(marima$y.analysed)
    rownames(Y) <- Names
    
    AR <- pol.mul( marima$ar.estimate, dif.poly,
    L  <- ( dim(marima$ar.estimate)[3] + dim(dif.poly)[3]-2 ) )
    MA <- marima$ma.estimate
    p  <- dim(AR)[3]
    q  <- dim(MA)[3]

    print(AR)
    print(MA)

    Rand <- rep(TRUE, d[1])
    Regr <- !Rand

    sig2 <- marima$resid.cov

        # print(sig2,digits=2)

    for (i in 1:d[1]) {
        if (sum(abs(AR[i, , ])) + sum(abs(MA[i, , ])) == 2) {
            Rand[i] <- FALSE
            sig2[i, ] <- 0
            sig2[, i] <- 0
        }
        if (sum(abs(AR[, i, ])) + sum(abs(MA[, i, ])) > 2) {
            Regr[i] <- TRUE
        }
    }

    L <- nstart + nstep
    if (L > d[2]) {
        cat("Input series length =", d[2], " but (nstart+nstep)=", L, "\n")
    }
        for (i in 1:d[1]) {
        # cat(Rand[i],Regr[i],'\n')
        if (!Rand[i] & Regr[i]) {
            cat("variable no.", i, "not random and regressor.\n")
            zy <- Y[i, (nstart:(L-1))]
            # zy[5]<-NaN
            if (is.na(sum(zy)) | is.nan(sum(zy))) {
        cat("Variable no.", i, ": illegal (NaN or NA) in future values. \n")
            }
            if (!is.na(sum(zy)) & !is.nan(sum(zy))) {
                cat("Variable no.", i, ": future values seem OK.\n")
            }
        }
    }
    # cat('Calling series : \n') print(Y[,(d[2]-19):d[2]])

    forecast <- arma.forecastingY(series = Y, ar.poly = AR,
       ma.poly = MA, nstart = nstart, nstep = nstep, Rand = Rand,
          Regr = Regr, marima = marima)

    pred.var <- forec.var(marima, nstep = nstep, dif.poly )$pred.var

    return(list(nstart = nstart, nstep = nstep,
        forecasts = forecast$estimates,
        residuals = forecast$residuals,
        pred.var = pred.var, averages = averages,
        mean.pattern = means))
}

arma.forecastingY <- function(series = NULL, ar.poly = NULL,
  ma.poly = NULL, nstart = NULL, nstep = NULL, Rand = NULL,
  Regr = NULL, marima = NULL) {
    "[" <- function(x, ...) .Primitive("[")(x, ..., drop = FALSE)
    d <- dim(series)
    if (d[1] > d[2]) {
        series <- t(series)
    }
    d <- dim(series)
    k <- d[1]
    N <- d[2]

    means    <- marima$mean.pattern
    averages <- marima$averages
    Constant <- marima$Constant
    cat("dimensions of time series", k, N, "\n")

    ar <- dim(ar.poly)
    ma <- dim(ma.poly)

    yseries <- matrix(c(series), nrow = d[1])

    yseries[is.na(yseries) | is.nan(yseries)] <- 0
    residuals <- yseries * 0
    forecasts <- residuals
    kvar <- ar[1]

    extra <- ar[3]
    su0 <- matrix(0, nrow=kvar, ncol=1)

    if( extra > 1 ){
        for (i in 1:extra){
            forecasts[, i] <- averages * means
            residuals[, i] <- yseries[, i] - forecasts[, i]           
        }
    }

    for (i in extra:N){
        suar <- su0
        if(ar[3] > 1) {
            for (j in 2:ar[3]) {
                suar <- suar + matrix(ar.poly[, , j], nrow = kvar) %*%
                    matrix(yseries[, (i+1-j)], ncol=1)
                if (i < 0)
                    cat("i, suar= ", i, "\n")
            }
        }

        suma <- su0
        if (ma[3] > 1) {
            for (j in 2:ma[3]) {
                if( i+1-j > 0 )  {
                suma <- suma + matrix(ma.poly[, , j], nrow = kvar) %*%
                    matrix(residuals[, (i + 1 - j)], ncol = 1)
                }
            }
        }

        est <- suma - suar + Constant
        forecasts[, i] <- est
        
        # cat('ar- & ma-polynomier \n')
        # print(ar.poly[,,2])
        # print(ma.poly[,,2])

        if (i <= nstart) {
            for (j in 1:k) {
                if (Rand[j]) {
                  residuals[j, i] <- yseries[j, i] - est[j]
                }
            }
        }
        if (i > nstart) {
            for (j in 1:k) {
                if (Rand[j]) {
                  yseries[j, i] <- est[j]
                  residuals[j, i] <- 0
                }
            }
        }
    }  # End i-loop here
    return(list(estimates = forecasts, residuals = residuals))
}
    
##' @title forec.var
##' 
##' @description Function for calculation of variances of nstep forecasts
##' using a marima type model.
##'
##' @param marima   =  marima object (cov.u and ar.estimates and
##'                    ma.estimates are used)
##' @param nstep    =  length of forecast
##'
##' @param dif.poly =  autoregressive representation of differencing
##'                    polynomial as constructed by the function
##'                    define.dif(...) when the time series is differenced
##'                    (if so) before being analysed by marima. 
##'
##' @return = pred.var   = variance-covariances for nstep forecasts
##' (an array with dimension (kvar,kvar,nstep).
##' 
##' @return = rand.shock = corresponding random shock representation
##' of the model used.
##'

forec.var <- function(marima, nstep = 1, dif.poly = NULL) {
    
    cat("Calculation forecasting variances.  \n")
    sig2 <- marima$resid.cov
    ar.poly <- marima$ar.estimates
    ma.poly <- marima$ma.estimates

    d <- dim(sig2)

    if(is.null(dif.poly)) {
        dif.poly <- diag(d[1])
    }
    
    dif.poly <- check.one(dif.poly)

    kvar = dim(sig2)[1]
    
    ar.poly<-pol.mul(ar.poly, dif.poly, L=nstep+1)

    if (is.null(ma.poly)) {
        ma.poly <- diag(kvar)
    }
    ar.poly <- check.one(ar.poly)
    ma.poly <- check.one(ma.poly)
    if (nstep < 1) {
        nstep <- 1
    }
    xsi.poly <- rand.shock(ar.poly, ma.poly, nstep)
    var <- xsi.poly
    # cat('nsteps=',1:nstep,'\n')
    for (i in 1:nstep) {
        # cat('i=',i,'\n')
        var[, , i] <- matrix(xsi.poly[, , i], nrow = kvar) %*%
            sig2 %*% t(matrix(xsi.poly[, , i], nrow = kvar))
        if (i > 1) {
            var[, , i] <- var[, , i] + var[, , (i - 1)]
        }
    }
    d <- dim(var)[3] - 1
    
    # cat('d=',d,'\n')
    # cat('xsi.poly=',1:nstep,'\n')
    # print(round(xsi.poly,4))
    # cat('var=','\n') print(round(var[,,1:d],4))
    
    return(list(pred.var = var[, , 1:d], rand.shock = xsi.poly[, , 1:d]))
}

##' @title arma.filter 
##' 
##' @description Filtering of (kvar-variate) time series with marima
##' type model.
##' 
##' Calculation of residuals and filtered values of timeseries using
##' a marima model.
##'
##' @param series matrix holding the kvar by n multivariate timeseries
##' (if (kvar > n) the series is transposed and a warning is given).
##'
##' @param ar.poly (kvar,kvar,p+1) array containing autoregressive matrix
##' polynomial model part. If the filtering is to be performed for
##' undifferenced data when the analysis (in marima) was done for differenced
##' data, the input array ar.poly should incorporate the ar-representation
##' of the differensing operation (using, for example:
##' ar.poly <- pol.mul(ar.estimate, dif.poly,
##' L = ( dim(ar.estimates)[3]+dim(dif.poly)[3])), where 'dif.poly'
##' was obtained when differencing the time series (using define.dif)
##' before analysing it with marima (giving the ar.estimate) . 
##'
##' @param ma.poly (kvar,kvar,q+1) array containing moving average matrix
##' polynomial model part.
##'
##' If a leading unity matrix is not included in the ar- and/or the ma-part
##' of the model this is automatically taken care of in the function
##' (in that case the dimensions of the model arrays used in arma.filter()
##'  are, respectively, (kvar,kvar,p+1) and (kvar,kvar,q+1)).
##'
##' @param means vector (length = kvar) indicating whether means are
##' subtracted or not (0/1). Default : means=1 saying that all means
##' are subtracted (equivalent to means = c(1,1,...,1)).
##'
##' @return estimates = estimated values for input series
##'
##' @return residuals = corresponding residuals
##'
##' @return averages = averages of variables in input series
##'
##' @return mean.pattern = pattern of means as used in filtering
##'
##' @examples
##'
##' library(marima)
##' data(austr)
##' series<-t(austr)[,1:90]
##' # Define marima model
##' Model5 <- define.model(kvar=7,ar=1,ma=1,rem.var=1,reg.var=6:7)
##'
##' # Estimate marima model
##' Marima5 <- marima(series,Model5$ar.pattern,Model5$ma.pattern,penalty=1)
##'
##' # Calculate residuals by filtering
##' Resid <- arma.filter(series,Marima5$ar.estimates,
##'      Marima5$ma.estimates)
##' # Compare residuals
##'
##' plot(Marima5$residuals[2,4:89], Resid$residuals[2,5:90],
##' xlab='marima residuals', ylab='arma.filter residuals')
##'
##' @export

arma.filter <- function(series = NULL, ar.poly = array(diag(kvar),
       dim = c(kvar, kvar, 1)), ma.poly = array(diag(kvar),
       dim = c(kvar, kvar, 1)), means = 1) {
    "[" <- function(x, ...) .Primitive("[")(x, ..., drop = FALSE)
    
    cat("arma.filter is being called \n")

    if (is.null(series)) {
        stop("No input series specified in input to arma.filter \n")
    }
    d <- dim(series)
    if (d[1] > d[2]) {
        cat("Warning: Input series is transposed to be a k x n series. \n")
        cat("Output (estimated values and residuals) \n")
        cat("will be organised the same way (k x n). \n")
        series <- t(series)
    }

    vmeans <- means
    # cat("vmeans,means=", vmeans, means, "\n")
    d <- dim(series)
    kvar <- d[1]

    averages <- rep(0, kvar)
    means <- rep(1, kvar)
    if (length(vmeans == 1)) {
        if (vmeans != 1) {
            means <- rep(0, kvar)
        }
    }
    for (i in 1:kvar) {
        averages[i] <- mean(series[i, ])
        if (means[i] == 1) {
            series[i, ] <- series[i, ] - averages[i]
        }
    }
    cat("indicators for means=", means, "\n")
    # cat(averages, "\n")

    ar.poly <- check.one(ar.poly)
    ma.poly <- check.one(ma.poly)
    
    m <- dim(ar.poly)
    ma <- dim(ma.poly)

    su0 <- matrix(0, nrow = kvar, ncol = 1)
    extra <- 2 * max(m[3], ma[3])

    extray <- matrix(0, nrow = kvar, ncol = extra)
    for (i in 1:kvar) {
        if (means[i] != 1) {
            extray[i, ] <- extray[i, ] + averages[i]
        }
    }

    # cat('extra=',extra,'\n') print(extray)

    yseries <- cbind(extray, series)
    estimates <- yseries * 0
    residuals <- estimates

    # cat('Series= \n') print(averages) print(round(yseries[,1:20],2))
    # print(round(short.form(ar.poly,leading=F),2))

    s <- dim(yseries)
    cat(" dim(yseries)", s, "\n")

    for (i in (extra/2 + 1):s[2]) {
        sur <- su0
        if (m[3] > 1) {
            for (j in 2:m[3]) {
                # cat('sur=',sur,'\n') cat('m,i,j=',m,i,j,'\n')
                # cat('dimensions=',dim(ar.poly[,,j]),dim(yseries[,(i+1-j)]),
                # '\n')
                # print(ar.poly[,,j]) print(yseries[,(i+1-j)])
                sur <- sur + matrix(ar.poly[, , j], nrow = kvar) %*%
                    matrix(yseries[, (i + 1 - j)], nrow = kvar)
            }
        }

        # cat(i,j,(i+1-j),sur,'\n') cat('y=',yseries,'\n')

        suma <- su0
        if (ma[3] > 1) {
            for (j in 2:ma[3]) {
                suma <- suma + matrix(ma.poly[, , j], nrow = kvar) %*%
                    matrix(residuals[, (i + 1 - j)], nrow = kvar)
            }
        }

        # cat(i,j,(i+1-j),suma,'\n')
        su <- -sur + suma
        estimates[, i] <- su
        residuals[, i] <- yseries[, i] - estimates[, i]
    }

    estimates <- estimates[, (extra + 1):(extra + d[2])]
    residuals <- residuals[, (extra + 1):(extra + d[2])]

    for (i in 1:kvar) {
        if (means[i] == 1) {
            estimates[i, ] <- estimates[i, ] + averages[i]
        }
    }

    # cat('dimensions',dim(yseries),dim(series),dim(residuals),'\n')

    return(list(estimates = estimates,
                residuals = residuals,
                averages = averages,
                mean.pattern = means))
}

    

