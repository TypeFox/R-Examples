.computeLength <- function(d) {
  res <- length(d)
  return(res)
}


.computeTrendStrength <- function(d) {
  # computing as a p-value of the zero slope in linear regression model
  # TODO: nebylo by lepsi vracet nejak logaritmus pval?)
  # TODO: prokladame primkou (linearni regrese) - co kdybychom zkusili prokladat taky necim jinym
  #       a definovat tak vice ruznych trendStrength-u?
  t <- 1:length(d)
  pval <- summary(lm(d ~ t))$coefficients[2, 4]
  return(1 - pval)
}


.computeSeasonStrength <- function(d) {
  # TODO: nebylo by lepsi vracet nejak logaritmus pval?)
  if (frequency(d) == 1) {
    # perioda = 12 mesicu, rocni casovka automaticky neni sezonni
    return(0)
  } else {
    time <- 1:length(d)
    num <- frequency(d) - 1
    vec <- rep(c(rep(0, num), 1), length(d))
    flags <- matrix(vec[1:(length(d) * num)], byrow=TRUE, ncol=num, nrow=length(d))
    form <- paste("d ~ time +", paste("flags[,", 1:num, "]", sep="", collapse=" + "))
    pval <- min(summary(lm(as.formula(form)))$coefficients[c(-1, -2), 4])
    return(1 - pval)
  }
}


.computeSkewness <- function(d) {
  # TODO: mozna by stalo za to  to transformovat nejak nelinearne
  s <- abs(skewness(d, type=1))
  return(s)
}


.computeKurtosis <- function(d) {
  # TODO: asi urcite by stalo za to  to transformovat nejak nelinearne
  k <- 3 + kurtosis(d, type=1)
  return(k)
}


.computeVarcoef <- function(d) {
  # TODO: uprava i pro zaporna data (tj. data s nulovym prumerem)
  return(sd(d) / mean(d))
}


.computeStationarity <- function(d) {
  # TODO: nebylo by lepsi vracet nejak logaritmus pval?)

  # on very short time series, this sometimes causes error
  # (see https://stackoverflow.com/questions/17282788/r-error-with-adf-test-in-time-series-lapply)
  # to recover that, we return 0.5
  pval <- try(suppressWarnings(adf.test(d)$p.value))
  if (inherits(pval, 'try-error')) {
      return(0.5)
  }
  return(1 - pval)
}


.computeFrequency <- function(d) {
  return(1 / frequency(d))
}





frbe <- function(d, h=10) {
    if (!is.ts(d)) {
        stop("'d' must be a time-series object")
    }

    if (!is.numeric(h) || length(h) != 1 || h < 1) {
        stop("'h' must contain a single positive integer value")
    }

    result <- list()
    result$data <- d
    
    result$forecasts <- data.frame(
                arima=as.numeric(forecast(auto.arima(d, stepwise=FALSE), h=h)$mean),
                expSmooth=as.numeric(forecast(ets(d), h=h)$mean), 
                randomWalk=as.numeric(rwf(d, drift=FALSE, h=h)$mean),
                theta=as.numeric(thetaf(d, h=h)$mean))


    result$features <- data.frame(length=.computeLength(d),
                                  trendStrength=.computeTrendStrength(d),
                                  seasonStrength=.computeSeasonStrength(d),
                                  skewness=.computeSkewness(d),
                                  kurtosis=.computeKurtosis(d),
                                  varcoef=.computeVarcoef(d),
                                  stationarity=.computeStationarity(d),
                                  frequency=.computeFrequency(d))
    f <- lcut3(result$features, context=.frbemodel$featuresContext)

    result$weights <- sapply(names(.frbemodel$model),
                             function(n) {
                                 ctx <- .frbemodel$weightContext[[n]]
                                 vals <- slices(ctx[1], ctx[3], 1000)
                                 parts <- lcut3(vals, name='weight', context=ctx)
                                 pbld(f, .frbemodel$model[[n]], parts, vals, type='global')
                             })
    result$weights <- result$weights[colnames(result$forecasts)]

    if (sum(result$weights) == 0) {
        result$weights <- rep(1, ncol(result$forecasts))
        names(result$weights) <- colnames(result$forecasts)
    }

    result$mean <- apply(result$forecasts, 1,
                         function(row) { 
                             sum(row * result$weights) / sum(result$weights) 
                         })

    class(result) <- c('frbe', class(result))
    return(result)
}
