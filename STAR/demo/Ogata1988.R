require(STAR)
originalpar <- par(ask = dev.interactive(orNone = TRUE))
oldpar <- par()
data(ShallowShocks)

## load ShallowShocks data
layout(matrix(1:3, nrow = 3))

## Reproduce Fig. 2
plot(ShallowShocks$Date,
     cumsum(ShallowShocks$energy.sqrt) / 10^13,
     type ="l",
     xlab = "",
     ylab = "",
     main = "Cumulative square root of energy")
plot(ShallowShocks$Date,
     cumsum(1+numeric(dim(ShallowShocks)[1])),
     type ="l",
     xlab = "",
     ylab = "",
     main = "Cumulative number of shocks")
plot(ShallowShocks$Date,
     ShallowShocks$magnitude,
     type = "h",
     ylim = c(5,9),
     xlab = "Time (days)",
     ylab = "",
     main = "Magnitude vs Occurrence time")

## Reproduce Fig. 3
layout(matrix(1))
magnitudeD.x <- sort(unique(ShallowShocks$magnitude))
magnitudeD.y <- as.vector(table(ShallowShocks$magnitude))
GutenbergRichter <- sapply(1:length(magnitudeD.y), function(idx) sum(magnitudeD.y[idx:length(magnitudeD.y)]))
plot(magnitudeD.x,
     GutenbergRichter,
     log = "y",
     xlab = "Magnitude",
     ylab = "Cumulative Number",
     main = "Check of Gutenberg-Richter Law")
points(magnitudeD.x, magnitudeD.y, pch = "*")

## Reproduce Fig. 4
## Get the sorted inter shock intervals
isi.sort <- sort(diff(ShallowShocks$Date))
survivor <- cumsum(1+numeric(length(isi.sort)))[length(isi.sort):1]
plot(isi.sort,
     survivor,
     pch = "+",
     log = "y",
     xlab = "Time Interval (days)",
     ylab = "Cumulative Number",
     main = "Empirical Log-Survivor Function"
     )

## Reproduce Fig. 5
vtShallow <- varianceTime(ShallowShocks$Date,,c(5,10,20,40,60,80,seq(100,500,by = 25))*10)
plot(vtShallow, style="Ogata")

## Reproduce Fig. 6
shalShocks <- lockedTrain(as.spikeTrain(ShallowShocks$Date),,c(0,500))
shalShocksH <- hist(shalShocks,5,plot=FALSE)
plot(shalShocksH,"Ogata",c(0.95,0.99),xlab="TIME LAG (DAYS)",ylab="NUMBER OF EVENTS PER DAY")

## Reproduce Fig. 7
myBinNb <- 101
myMidPoints <- seq(from = 1, to = 6, length.out=myBinNb)
myMidPoints <- 10^myMidPoints/200
myBreaks <- c(0,myMidPoints[-myBinNb] + diff(myMidPoints) / 2)
shalShocksH2 <- hist(shalShocks,breaks=myBreaks,plot=FALSE)
yy <- abs(shalShocksH2$density - shalShocksH2$refFreq)
plot(shalShocksH2$mids[shalShocksH2$density>0],
     yy[shalShocksH2$density>0],
     pch = 1,
     xlim = c(0.001,10000),
     log = "xy",
     xlab = "TIME LAG (DAYS)",
     ylab = "NUMBER OF EVENTS PER DAY"
     )

##################################################################
## Define a couple of functions for parametric model fit #########
##################################################################

makeCI.Omori <- function(data = ShallowShocks,
                         modelType = "Epidemic",
                         mu.initial = 0.00536,
                         K.initial = 0.017284,
                         c.initial = 0.01959,
                         p.initial = 1.0,
                         beta.initial = 1.61385
                         ) {

  ## check that the first argument is a data frame
  if (!identical(class(data), "data.frame")) stop("data should be a data frame.")
  ## check that data has a variable called Date
  if (!"Date" %in% names(data)) {
    stop("data should contained a Date variable.")
  } else {
    Date <- data$Date
  }
  ## if magnitudeEffect is true, check that data contains a magnitude variable 
  if (beta.initial != 0) {
    if (! "magnitude" %in% names(data)) {
      stop("data should contain a magnitude variable.")
    } else {
      ## create a magnitude variable
      ## We will assume that the minimal magnitude value corresponds
      ## to the cut off used. We the offset every magnitude with
      ## respect to this cut off.
      magnitude <- data$magnitude - min(data$magnitude)
    }
  }
  ## if modelType is Trigger, check that data contains a type variable
  ## if it does, create a type variable and set elements with "foreshock"
  ## value to "main"
  if (modelType == "Trigger") {
    if (! "type" %in% names(data)) {
      stop("data should contain a type variable.")
    } else {
      type <- data$type
      type[type == "foreshock"] <- "main"
    }
  }

  ## We do not need data anymore beyond this point so we get rid of it
  rm(data)

  Eq.14 <- function(t, K, c, p) K / (t + c)^p

  Eq.14.sum <- function(t, K, c, p) {
    if (p != 1) return(K / (1 - p) * ((t + c)^(1-p) - c^(1-p)))
    else return(K * log((t+c)/c))
  }

  if (beta.initial != 0) { ## magnitude effects are taken into account

    CI <- function(t, mu, K, c, p, beta) {
        if (missing(mu)) mu <- mu.initial
        if (missing(K)) K <- K.initial
        if (missing(c)) c <- c.initial
        if (missing(p)) p <- p.initial
        if (missing(beta)) beta <- beta.initial
        result <- mu
        if (modelType == "Trigger") {
          goodTimes <- t - Date[type == "main" & Date < t]
          goodEffects <- exp(beta * magnitude[type == "main" & Date < t])
          result <- result + ifelse(length(goodTimes) > 0, sum(goodEffects * sapply(goodTimes, Eq.14, K = K, c = c, p = p)), 0)
        } else {
          goodTimes <- t - Date[Date < t]
          goodEffects <- exp(beta * magnitude[Date < t])
          result <- result + ifelse(length(goodTimes) > 0, sum(goodEffects * sapply(goodTimes, Eq.14, K = K, c = c, p = p)), 0)
        } ## End of the conditional on modelType == "Trigger"
        return(result)
      } ## End of CI definition

    CIsum <- function(t, mu, K, c, p, beta) {
        if (missing(mu)) mu <- mu.initial
        if (missing(K)) K <- K.initial
        if (missing(c)) c <- c.initial
        if (missing(p)) p <- p.initial
        if (missing(beta)) beta <- beta.initial
        result <- mu * t
        if (modelType == "Trigger") {
          goodTimes <- t - Date[type == "main" & Date < t]
          goodEffects <- exp(beta * magnitude[type == "main" & Date < t])
          result <- result + ifelse(length(goodTimes) > 0, sum(goodEffects * sapply(goodTimes, Eq.14.sum, K = K, c = c, p = p)), 0)
        } else {
          goodTimes <- t - Date[Date < t]
          goodEffects <- exp(beta * magnitude[Date < t])
          result <- result + ifelse(length(goodTimes) > 0, sum(goodEffects * sapply(goodTimes, Eq.14.sum, K = K, c = c, p = p)), 0)
        } ## End of the conditional on modelType == "Trigger"
        return(result)
      } ## End of CIsum definition

    } else {

      CI <- function(t, mu, K, c, p) {
        if (missing(mu)) mu <- mu.initial
        if (missing(K)) K <- K.initial
        if (missing(c)) c <- c.initial
        if (missing(p)) p <- p.initial
        result <- mu
        if (modelType == "Trigger") {
          goodTimes <- t - Date[type == "main" & Date < t]
          result <- result + sum(sapply(goodTimes, Eq.14, K = K, c = c, p = p))
        } else {
          goodTimes <- t - Date[Date < t]
          result <- result + sum(sapply(goodTimes, Eq.14, K = K, c = c, p = p))
        } ## End of the conditional on modelType == "Trigger"
        return(result)
      } ## End of CI definition

      CI <- function(t, mu, K, c, p) {
        if (missing(mu)) mu <- mu.initial
        if (missing(K)) K <- K.initial
        if (missing(c)) c <- c.initial
        if (missing(p)) p <- p.initial
        result <- mu * t
        if (modelType == "Trigger") {
          goodTimes <- t - Date[type == "main" & Date < t]
          result <- result + sum(sapply(goodTimes, Eq.14.sum, K = K, c = c, p = p))
        } else {
          goodTimes <- t - Date[Date < t]
          result <- result + sum(sapply(goodTimes, Eq.14.sum, K = K, c = c, p = p))
        } ## End of the conditional on modelType == "Trigger"
        return(result)
      } ## End of CIsum definition
      
    } ## End of conditional on magnitudeEffect

  return(list(call = match.call(),
              CI = CI,
              CIsum = CIsum,
              modelType = modelType
              )
         )
}


makeMinusLogLik.Omori <- function(data = ShallowShocks,
                                  modelType = "Epidemic",
                                  withBeta = TRUE,
                                  withP = TRUE,
                                  mu.initial = 0.00536,
                                  K.initial = 0.017284,
                                  c.initial = 0.01959,
                                  p.initial = 1.0,
                                  beta.initial = 1.61385,
                                  observationStart = 0,
                                  observationEnd = as.numeric(as.Date("1980-12-31") - as.Date("1885-1-1"))
                                  ) {

  ## check that the first argument is a data frame
  if (!identical(class(data), "data.frame")) stop("data should be a data frame.")
  ## check that data has a variable called Date
  if (!"Date" %in% names(data)) {
    stop("data should contained a Date variable.")
  } else {
    Date <- data$Date
  }
  ## if magnitudeEffect is true, check that data contains a magnitude variable 
  if (withBeta) {
    if (! "magnitude" %in% names(data)) {
      stop("data should contain a magnitude variable.")
    } else {
      ## create a magnitude variable
      ## We will assume that the minimal magnitude value corresponds
      ## to the cut off used. We the offset every magnitude with
      ## respect to this cut off.
      magnitude <- data$magnitude - min(data$magnitude)
    }
  }
  ## if modelType is Trigger, check that data contains a type variable
  ## if it does, create a type variable and set elements with "foreshock"
  ## value to "main"
  if (modelType == "Trigger") {
    if (! "type" %in% names(data)) {
      stop("data should contain a type variable.")
    } else {
      type <- data$type
      type[type == "foreshock"] <- "main"
    }
  }

  ## We do not need data anymore beyond this point so we get rid of it
  rm(data)

  Eq.14 <- function(t, K, c, p) K / (t + c)^p

  Eq.14.sum <- function(t, K, c, p) {
    if (p != 1) return(K / (1 - p) * ((t + c)^(1-p) - c^(1-p)))
    else return(K * log((t+c)/c))
  }

  CI <- function(t, mu, K, c, p, beta) {
    if (missing(mu)) mu <- mu.initial
    if (missing(K)) K <- K.initial
    if (missing(c)) c <- c.initial
    if (missing(p)) p <- ifelse(withP, p.initial, 1.0)
    if (missing(beta)) beta <- ifelse(withBeta, beta.initial, 0.0)
    if (modelType == "Trigger") {
      goodTimes <- t - Date[type == "main" & Date < t]
      if (withBeta) goodEffects <- exp(beta * magnitude[type == "main" & Date < t])
      else goodEffects <- 1
      result <- ifelse(length(goodTimes) > 0, sum(goodEffects * sapply(goodTimes, Eq.14, K = K, c = c, p = p)), 0)
      return(result)
    } else {
      goodTimes <- t - Date[Date < t]
      if (withBeta) goodEffects <- exp(beta * magnitude[Date < t])
      else goodEffects <- 1
      result <- mu + ifelse(length(goodTimes) > 0, sum(goodEffects * sapply(goodTimes, Eq.14, K = K, c = c, p = p)), 0)
      return(result)
    } ## End of the conditional on modelType == "Trigger"
    
  } ## End of CI definition
  
  CIsum <- function(t.start, t.end, mu, K, c, p, beta) {
    if (missing(t.start)) t.start <- observationStart
    if (missing(t.end)) t.end <- observationEnd
    if (missing(mu)) mu <- mu.initial
    if (missing(K)) K <- K.initial
    if (missing(c)) c <- c.initial
    if (missing(p)) p <- ifelse(withP, p.initial, 1.0)
    if (missing(beta)) beta <- ifelse(withBeta, beta.initial, 0.0)
    result <- mu * (t.end - t.start)
    if (modelType == "Trigger") {
      goodTimes <- t.end - Date[type == "main" & Date < t.end & Date >= t.start]
      if (withBeta) goodEffects <- exp(beta * magnitude[type == "main" & Date < t.end & Date >= t.start])
      else goodEffects <- 1
      goodEffects <- ifelse(withBeta, exp(beta * magnitude[type == "main" & Date < t.end & Date >= t.start]), 1)
      result <- result + ifelse(length(goodTimes) > 0, sum(goodEffects * sapply(goodTimes, Eq.14.sum, K = K, c = c, p = p)), 0)
    } else {
      goodTimes <- t.end - Date[Date < t.end & Date >= t.start]
      if (withBeta) goodEffects <- exp(beta * magnitude[Date < t.end & Date >= t.start])
      else goodEffects <- 1
      result <- result + ifelse(length(goodTimes) > 0, sum(goodEffects * sapply(goodTimes, Eq.14.sum, K = K, c = c, p = p)), 0)
    } ## End of the conditional on modelType == "Trigger"
    return(result)
  } ## End of CIsum definition
  
  minusLogLik <- function(theta) {
    paraNames <- c("mu", "K", "c", "p", "beta")
    if (!all(names(theta) %in% paraNames)) stop("Some parameters names are not right.")
    if (!withP) theta["p"] <- 1.0
    if (!withBeta) theta["beta"] <- 0.0
    ## All parameters should be positive or null
    if (any(theta < 0)) return(Inf)
    paraList <- as.list(theta)
    ## Get first the contribution of the integral of the intensity
    noEventContrib <- do.call(CIsum, c(list(t.start = observationStart, t.end = observationEnd), paraList))
    ## Get the contibutions of the events
    if (modelType != "Trigger") { 
      eventContrib <- sum( sapply(Date, function(d) -log(do.call(CI, c(list(t = d), paraList))) ) )
    } else {
      eventContrib.a <- sum( sapply(Date[type != "main"], function(d) -log(do.call(CI, c(list(t = d), paraList))) ) )
      eventContrib.c <- - sum(type == "main") * log(theta["mu"]) 
      eventContrib <- eventContrib.a + eventContrib.c
    }
    return( noEventContrib + eventContrib )
  }

  theta.initial <- c(mu = mu.initial,
                     K = K.initial,
                     c = c.initial)

  if (withP) theta.initial["p"] <- p.initial
  if (withBeta) theta.initial["beta"] <- beta.initial
  
  return(list(call = match.call(),
              CI = CI,
              CIsum = CIsum,
              minusLogLik = minusLogLik,
              theta.initial = theta.initial,
              modelType = modelType,
              withP = withP,
              withBeta = withBeta,
              observationStart = observationStart,
              observationEnd = observationEnd
              )
         )
}

############################################################
## Get the bottom right part of Table 2
############################################################

## Define model with free beta and p fixed at 1
omoriBetaNoPepidemic <- makeMinusLogLik.Omori(withP = FALSE)
## The following commands takes 208 s on pdp8, 276 s on a laptop PIV 3 GHz
omoriBetaNoPepidemic.MLE1 <- optim(omoriBetaNoPepidemic$theta.initial,
                                   omoriBetaNoPepidemic$minusLogLik,
                                   control = list(REPORT=1, trace = 2),
                                   hessian = TRUE,
                                   method="BFGS")
## Check the etsimates
rbind(omoriBetaNoPepidemic.MLE1$par, sqrt(diag(solve(omoriBetaNoPepidemic.MLE1$hessian))))

## Define model with free beta and p
omoriBetaPepidemic <- makeMinusLogLik.Omori()
## The following commands takes 415 s on pdp8, 588 s on a laptop PIV 3 GHz
omoriBetaPepidemic.MLE1 <- optim(omoriBetaPepidemic$theta.initial,
                                 omoriBetaPepidemic$minusLogLik,
                                 control = list(REPORT=1, trace = 2),
                                 method="BFGS")

## Define model with beta fixed at 0 and p fixed at 1
omoriNoBetaNoPepidemic <- makeMinusLogLik.Omori(withBeta = FALSE, withP = FALSE)
## The following commands takes 132 s on pdp8, 195 s on a laptop PIV 3 GHz
omoriNoBetaNoPepidemic.MLE1 <- optim(omoriBetaNoPepidemic.MLE1$par[c("mu","K","c")],
                                     omoriNoBetaNoPepidemic$minusLogLik,
                                     control = list(REPORT=1, trace = 2),
                                     method="BFGS")

## Define model with beta fixed at 0 and p free
omoriNoBetaPepidemic <- makeMinusLogLik.Omori(withBeta = FALSE, withP = TRUE)
## The following commands takes 143 s on pdp8, 186 s on a laptop PIV 3 GHz
omoriNoBetaPepidemic.MLE1 <- optim(c(omoriNoBetaNoPepidemic.MLE1$par, p = 1.0),
                                   omoriNoBetaPepidemic$minusLogLik,
                                   control = list(REPORT=1, trace = 2),
                                   method="BFGS")

## show minus log lik, number of parameters and AIC
myMinusLogLik <- c(omoriNoBetaNoPepidemic.MLE1$value,
                   omoriNoBetaPepidemic.MLE1$value,
                   omoriBetaNoPepidemic.MLE1$value,
                   omoriBetaPepidemic.MLE1$value)
myNbPar <- as.integer(c(length(omoriNoBetaNoPepidemic.MLE1$par),
                        length(omoriNoBetaPepidemic.MLE1$par),
                        length(omoriBetaNoPepidemic.MLE1$par),
                        length(omoriBetaPepidemic.MLE1$par))
                      )
myAIC <- 2 * myMinusLogLik + 2 * myNbPar
mySummary <- rbind(myMinusLogLik,
                   myNbPar,
                   myAIC)
rownames(mySummary) <- c("-log L(theta)",
                         "Number of Parameters",
                         "AIC")
colnames(mySummary) <- c("beta = 0, p = 1",
                         "beta = 0, p free",
                         "beta free, p = 1",
                         "beta free, p free")
mySummary
############################################################
## Bottom right part of Table 2 done
############################################################                   

## Replicate Fig. 8 of Ogata 1988
## We build the Conditional Intensity with a resolution of 1 day
myTime <- seq(from=floor(min(ShallowShocks$Date)),
              to=ceiling(max(ShallowShocks$Date)),
              by=10
              )
## The next computation takes 173 s on PIV 3GHz
bestCI <- sapply(myTime,
                 omoriBetaNoPepidemic$CI,
                 mu = omoriBetaNoPepidemic.MLE1$par["mu"],
                 K = omoriBetaNoPepidemic.MLE1$par["K"],
                 c = omoriBetaNoPepidemic.MLE1$par["c"],
                 beta = omoriBetaNoPepidemic.MLE1$par["beta"],
                 p = 1.0)

plot(myTime,
     bestCI,
     type = "l",
     log = "y",
     xlab = "Time (days)",
     ylab = "Shocks / Day",
     main = "Estimated Conditional Intensity with Best Model"
     )
rug(ShallowShocks$Date)

## Replicate Fig. 9
bestLambda <- sapply(ShallowShocks$Date,
                     function(t) {
                       omoriBetaNoPepidemic$CIsum(t.start = 0,
                                                  t.end = t,
                                                  mu = omoriBetaNoPepidemic.MLE1$par["mu"],
                                                  K = omoriBetaNoPepidemic.MLE1$par["K"],
                                                  c = omoriBetaNoPepidemic.MLE1$par["c"],
                                                  beta = omoriBetaNoPepidemic.MLE1$par["beta"],
                                                  p = 1.0)
                     }
                     )

bestLambda <- mkCPSP(bestLambda)
bestLambda.summary <- summary(bestLambda)
plot(bestLambda.summary,which=1)


## Replicate Fig. 10
plot(bestLambda.summary,which=2)

## Replicate Fig. 11
plot(bestLambda.summary,which=4)

## Replicate Fig. 12
plot(bestLambda.summary,which=3)

## Replicate Fig. 13
plot(bestLambda.summary,which=5)
## the same but better (with more points)
plot(varianceTime(bestLambda$ppspFct(),,seq(2.5,70,2.5)),style="Ogata")

par(originalpar)
