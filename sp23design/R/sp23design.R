###################################################
### code chunk number 19: Analyze-Design
###################################################
analyzeSP23Design <- function(sp23Design, trialHistory = NULL, data = NULL,
                              col=c("red", "red", "brown", "brown"),
                              lty=c(1,2,1,2)) {
  result <- list(responseSummary=NA, designSummary= NA)
  trialParameters <- sp23Design$trialParameters
  numLooks <- length(trialParameters$interimLookTime)
  trialExitTime <- function(historyDF) {
    j <- which(historyDF$acceptH0R > 0)
    if (length(j) <= 0) { ## we did not accept H_0^R
      j <- which(historyDF$rejectH0S > 0)
      if (length(j) <= 0) { ## H_0^S and hence null was NOT rejected
        j <- which(historyDF$acceptH0S > 0) ## H_0^S was accepted somewhere
      }
    }
    trialParameters$interimLookTime[j] ## there should only be 1 such entry
  }
  trialExitLook <- function(historyDF) {
    j <- which(historyDF$acceptH0R > 0)
    if (length(j) <= 0) { ## we did not accept H_0^R
      j <- which(historyDF$rejectH0S > 0)
      if (length(j) <= 0) { ## H_0^S and hence null was NOT rejected
        j <- which(historyDF$acceptH0S > 0) ## H_0^S was accepted somewhere
      }
    }
    j
  }

  if (is.data.frame(trialHistory)) {
    if (is.null(data)) stop("Expect clinical trial data as parameter!")
    exitTime <- trialExitTime(trialHistory)
    interimData <- generateInterimData(data,
                                       exitTime,
                                       trialParameters$adminCensoringTime)
    m <- nrow(interimData)
    arm <- factor(ifelse(interimData$treatmentIndicator == 1, "treatment", "control"))
    responder <- factor(ifelse(interimData$responseIndicator == 1, "responder", "non-responder"))
    interimData$arm <- arm
    interimData$responder <- responder
    sFit <- survfit(Surv(timeToEvent, delta) ~ arm + responder, data = interimData)
    plot(sFit, col=col, lty=lty, lwd=2, xlab="Time to Event (yrs)", ylab="Survival",
         xlim=c(0, trialParameters$adminCensoringTime))
    legend <- c(paste(levels(arm)[1], "-", levels(responder)), paste(levels(arm)[2], "-", levels(responder)))
    legend(0.5, 0.2, legend, text.col=col, pt.lwd=2, lty=lty, col=col)
    result$responseSummary<- computeResponseSummary(interimData)
  } else {
    ##
    ## Report on the statistics
    ##
    numberOfTimesH0RIsRejectedAtFirstLook <- sum(sapply(trialHistory, function(x) x$rejectH0R[1]))
    numberOfTimesH0RIsRejected <- sum(sapply(trialHistory, function(x) (sum(x$rejectH0R) > 0)))
    numberOfTimesStoppedForFutility <- sum(sapply(trialHistory, function(x) sum(x$acceptH0R | x$acceptH0S)))
    numberOfTimesH0SIsRejected <- sum(sapply(trialHistory, function(x) (sum(x$rejectH0S) > 0)))
    numberOfTimesH0SIsAccepted <- sum(sapply(trialHistory, function(x) (sum(x$acceptH0S) > 0)))
    numberOfTimesFutilityDecidedAtLastLook <- sum(sapply(trialHistory, function(x) x$acceptH0S[numLooks]))
    numberOfTimesTrialEndedAtLook <- sapply(1:numLooks,
                                            function(j) sum(sapply(trialHistory, function(x) x$acceptH0R[j] || x$acceptH0S[j] || x$rejectH0S[j])))
    names(numberOfTimesTrialEndedAtLook) <- as.character(trialParameters$interimLookTime)
    trialExitTimes <- sapply(trialHistory, trialExitTime)
    avgExitTime <- mean(trialExitTimes)
    designSummary <- list(numberOfTimesH0RIsRejectedAtFirstLook = numberOfTimesH0RIsRejectedAtFirstLook,
                          numberOfTimesH0RIsRejected = numberOfTimesH0RIsRejected,
                          numberOfTimesStoppedForFutility = numberOfTimesStoppedForFutility,
                          numberOfTimesH0SIsAccepted = numberOfTimesH0SIsAccepted,
                          numberOfTimesH0SIsRejected = numberOfTimesH0SIsRejected,
                          numberOfTimesFutilityDecidedAtLastLook = numberOfTimesFutilityDecidedAtLastLook,
                          numberOfTimesTrialEndedAtLook = numberOfTimesTrialEndedAtLook,
                          avgExitTime=avgExitTime)
    result$designSummary <- designSummary
  }
  result
}


###################################################
### code chunk number 12: Boundary-Function-for-b
###################################################
mHP.b <- function(mu=c(0,0), v=c(1,2), alpha=0.05, eps=1/2, side=c("one", "two")) {
  side <- match.arg(side)
  vsq = sqrt(v)
  v0 = diag(v)
  k=length(v)
  if (k > 1) {
  for (i in 1:(k-1))
    for (j in (i+1):k)
      v0[i,j] = v0[j,i] = v[i];
  }

  alphaEps <- alpha * eps
  ## Calculate the prob of stopping at interim looks
  crossingProbDiff <- function(b=3.0) {
    if (length(v) == 1) {
      if (side == "one")
        pnorm(-b * vsq, mean=mu, sd=vsq) - alphaEps
      else
        2 * pnorm(-b * vsq, mean=mu, sd=vsq) - alphaEps
    } else {
      bound <- b * vsq
      if (side == "one")
        1 - pmvnorm(lower = -Inf, upper=bound, mean=mu, sigma=v0, algorithm=Miwa()) - alphaEps
      else
        1 - pmvnorm(lower = -bound, upper=bound, mean=mu, sigma=v0, algorithm=Miwa()) - alphaEps
    }
  }

  if (k==1) {
    if (side=="one")
      list(alpha=alpha, eps=eps, mu=mu, v=v, b=qnorm(1-alphaEps), alpha.spent=alphaEps)
    else if(side=="two")
      list(alpha=alpha, eps=eps, mu=mu, v=v, b=qnorm(1-alphaEps/2), alpha.spent=alphaEps)
  } else {
    w <- uniroot(f=crossingProbDiff,
                 lower=qnorm(1 - alphaEps/2),
                 upper=qnorm(1 - alphaEps / length(mu)/2))
    ##   cat("b=",w$root,"\n")
    list(alpha=alpha, eps=eps, mu=mu, v=v, b=w$root, alpha.spent=w$f.root + alphaEps)
  }
}


###################################################
### code chunk number 13: Boundary-Function-for-c
###################################################
mHP.c <- function(mu=c(0, 0, 0), v=c(1, 2, 3), b=3.0, alpha=0.05, eps=1/2, side=c("one", "two")) {
  side = match.arg(side)
  vsq = sqrt(v)
  bVsq = b * vsq
  k=length(v)
  v0=diag(v)
  for(i in 1:(k-1))
    for(j in (i+1):k)
      v0[i,j]=v0[j,i]=v[i];

  crossingProbDiff <- function(c1=3.0) {
    bound <- c(bVsq[1:(k-1)], c1*sqrt(v[k]))
    if (side == "one")
      1-pmvnorm(lower=-Inf, upper=bound, mean=mu, sigma=v0, algorithm=Miwa()) - alpha
    else
      1-pmvnorm(lower=-bound, upper=bound, mean=mu, sigma=v0, algorithm=Miwa()) - alpha
  }

  w <- uniroot(f=crossingProbDiff,
               lower=qnorm(1 - alpha),
               upper=b)
  list(alpha=alpha, eps=eps, mu=mu, v=v, b=b, c=w$root,
       type1=w$f.root+alpha)
}


###################################################
### code chunk number 5: ComputeDGivenXi
###################################################
computeDGivenXi <- function(piVec, xiVec) {
  pi0 <- piVec[1]; pi1 <- piVec[2]
  (pi0 * xiVec[1] + 1 - pi0) - (pi1 * xiVec[1] * xiVec[3] + 1 - pi1) * xiVec[2]
}


###################################################
### code chunk number 10: ComputeGammaSubT
###################################################
computeGammaSubT <- function(thetaHat, pi, interimData) {
  hes = hessian(thetaHat, pi, interimData)
  hes[3, 3] - t(hes[3, 1:2]) %*% solve(hes[1:2, 1:2]) %*% hes[1:2, 3]
}


###################################################
### code chunk number 3: ComputeResponseSummary-Function
###################################################
computeResponseSummary <- function(interimData) {
  m = nrow(interimData)
  interimNumberOfEvents = sum(interimData$delta)
  control <- which(interimData$treatmentIndicator == 0)  ## indices of those on control
  treatment <- which(interimData$treatmentIndicator == 1)  ## indices of those on treatment
  m0 <- length(control)  ## number on control
  m1 <- length(treatment)  ## number on treatment
  y0 <- sum(interimData$responseIndicator[control]) ## number of responses in control
  y1 <- sum(interimData$responseIndicator[treatment]) ## number of responses in treatment
  c(m = m, m0 = m0, y0 = y0, m1 = m1, y1 = y1, numberOfTotalResponses = y0 + y1,
    controlRespProp = y0 / m0, treatmentRespProp = y1 / m1, pooledProp = (y0 + y1) / (m0 + m1))
}


###################################################
### code chunk number 17: Design-Execution-Function
###################################################
executeSP23Design <- function(sp23DesignObject, data, currentCalendarTime) {
  ##Find the largest interim look time smaller than calendarTime
  interimLookTime <- sp23DesignObject$trialParameters$interimLookTime
  l <- which(interimLookTime <= currentCalendarTime)
  llen <- length(l)
  if (llen > 0) {
    k <- l[llen] # the last index
    trialParameters <- sp23DesignObject$trialParameters
    glrBoundary <- sp23DesignObject$glrBoundary
    interimLookHistoryDF = sp23DesignObject$interimLookHistoryDF
    interimData <- generateInterimData(data,
                                       trialParameters$interimLookTime[k],
                                       trialParameters$adminCensoringTime)
    ## this ordering is important for the computation of the
    ## likelihood functions. The ordering is by time to event.
    interimData <- interimData[order(interimData$timeToEvent), ]

    performInterimLook(k, sp23DesignObject$trueParameters, trialParameters, glrBoundary,
                       interimData, interimLookHistoryDF,
                       argRejectH0R = ifelse(k == 1, FALSE, interimLookHistoryDF$rejectH0R[k-1]))
  } else {
    stop(cat("Trial design does not permit an interim look at time", currentCalendarTime, "!\n"))
  }
}


###################################################
### code chunk number 15: Design-Generator-function
###################################################
generateSP23Design <- function(trueParameters, trialParameters) {
  interimLookTime <- trialParameters$interimLookTime
  ##
  ## Some sanity checks
  ##
  ## Interim look times must be ordered
  if (sum(abs(sort(interimLookTime) - interimLookTime)) > 0) {
    stop(cat("Trial interim Look times", interimLookTime, "must be ordered\n"))
  }
  numLooks <- length(interimLookTime)
  N = sum(trialParameters$numberRecruitedEachYear) ## total no of subjects

  ## Compute the boundaries for efficacy for response
  test1 = mHP.b(mu = numeric(numLooks - 2),
    v = cumsum(trialParameters$numberRecruitedEachYear[1:(numLooks - 2)]),  ## why 1:3???
    alpha = trialParameters$type1ErrorForResponse,
    eps = trialParameters$epsType1,
    side = trialParameters$glrBoundarySidedness)

  test2 = mHP.c(mu = numeric(numLooks - 1),
    v = cumsum(trialParameters$numberRecruitedEachYear),  ## why 1:3???,
    b = test1$b,
    alpha = trialParameters$type1ErrorForResponse,
    eps = trialParameters$epsType1,
    side = trialParameters$glrBoundarySidedness)

  b.resp = c(rep(test1$b, numLooks - 2), test2$c, Inf) ## efficacy boundary

  if (is.na(trialParameters$type2ErrorForResponse)) {
    b.tilde.resp <- rep(Inf, numLooks)
  } else {
    test3 = mHP.b(mu = numeric(numLooks - 2),
      v = cumsum(trialParameters$numberRecruitedEachYear[1:(numLooks - 2)]),
      alpha = trialParameters$type2ErrorForResponse,
      eps = trialParameters$epsType2,
      side = trialParameters$glrBoundarySidedness)
    b.tilde.resp = c(rep(test3$b, numLooks - 2), Inf, Inf) ## futility boundary
  }

  test1 = mHP.b(mu = numeric(numLooks - 1),
    v = 1:(numLooks - 1),
    alpha = trialParameters$type1Error,
    eps = trialParameters$epsType1,
    side = trialParameters$glrBoundarySidedness)
  b.metas = c(rep(test1$b, numLooks - 1), Inf)  ## efficacy boundary

  if (is.na(trialParameters$type2Error)) {
    b.tilde.metas <- rep(Inf, numLooks)
  } else {
    test2 = mHP.b(mu = rep(0, 4),
      v = 1:(numLooks - 1),
      alpha = trialParameters$type2Error,
      eps = trialParameters$epsType2,
      side = trialParameters$glrBoundarySidedness)
    b.tilde.metas = c(rep(test2$b, numLooks - 1), Inf)  ## futility boundary
  }

##  cat("b.resp=",b.resp,"\n")
##  cat("b.tilde.resp=",b.tilde.resp,"\n")
##  cat("b.metas=",b.metas,"\n")
##  cat("b.tilde.metas=",b.tilde.metas,"\n\n")
  glrBoundary <- matrix(c(b.resp, b.tilde.resp, b.metas, b.tilde.metas), ncol=4)
  colnames(glrBoundary) <- c("RespEfficacy", "RespFutility", "SurvEfficacy", "SurvFutility")


  naVector <- rep(NA, numLooks) # a qty we reuse
  logicalVector <- logical(numLooks) # another qty we reuse.
  interimLookHistoryDF <- data.frame(m0 = naVector,
                                     m1 = naVector,
                                     y0 = naVector,
                                     y1 = naVector,
                                     pi0Hat = naVector,
                                     pi1Hat = naVector,
                                     pi0HatH0 = naVector,
                                     pi1HatH0 = naVector,
                                     pi0HatH1 = naVector,
                                     pi1HatH1 = naVector,
                                     glrRespH0 = naVector,
                                     glrRespH1 = naVector,
                                     hazard = naVector,
                                     alphaHat = naVector,
                                     betaHat = naVector,
                                     gammaHat = naVector,
                                     alphaHatH0 = naVector,
                                     betaHatH0 = naVector,
                                     gammaHatH0 = naVector,
                                     alphaHatH1 = naVector,
                                     betaHatH1 = naVector,
                                     gammaHatH1 = naVector,
                                     glrSurvH0 = naVector,
                                     glrSurvH1 = naVector,
                                     v = naVector,
                                     usedV = naVector,
                                     rejectH0R = logicalVector,
                                     acceptH0R = logicalVector,
                                     rejectH0S = logicalVector,
                                     acceptH0S = logicalVector
                                     )
  list(trueParameters = trueParameters, trialParameters=trialParameters,
       glrBoundary=glrBoundary, interimLookHistoryDF=interimLookHistoryDF)
}


###################################################
### code chunk number 18: Explore-Design-Function
###################################################
exploreSP23Design <- function(sp23Design, numberOfSimulations = 25,
                              rngSeed = 12345, showProgress=TRUE) {
  set.seed(rngSeed)
  trialHistory <-  vector(mode="list", length=numberOfSimulations)
  if (showProgress) pb <- txtProgressBar(min = 0, max = numberOfSimulations, style = 3)
  for (i in 1:numberOfSimulations) {
    sp23DesignObj <- resetSP23Design(sp23Design) ## reset statistics
    trialParameters <- sp23DesignObj$trialParameters
    trueParameters <- sp23DesignObj$trueParameters
    d <- generateClinicalTrialData(nRec = trialParameters$numberRecruitedEachYear,
                                   nFUp = trialParameters$followupTime,
                                   pi0 = trueParameters$p0,
                                   pi1 = trueParameters$p1,
                                   theta = trueParameters$theta,
                                   lambda0 = trueParameters$baselineLambda)
    interimLookTime <- sp23DesignObj$trialParameters$interimLookTime
    numLooks <- length(interimLookTime)
    for (k in 1:numLooks) {
      result <- executeSP23Design(sp23DesignObj, d, interimLookTime[k])
      ## We gather the data for the interim history for the next look
      for (x in names(sp23DesignObj$interimLookHistoryDF)) {
        sp23DesignObj$interimLookHistoryDF[k, x] <- result[x]
      }
      ## Need to also store the last beta if calculated, but easier to store it any way
      if (k == numLooks) {
        sp23DesignObj$glrBoundary[k, "SurvEfficacy"] <- result["b.metas.Last"]
        sp23DesignObj$glrBoundary[k, "SurvFutility"] <- result["b.metas.Last"]
      }
      ##print(result)
      if (result["acceptH0R"] || result["acceptH0S"] ||
          (result["rejectH0R"] && result["rejectH0S"]) ) {
        ## we accepted H_0^R, or accepted H_0^S or rejected the composite null of (H_0^R & H_0^S)
        ## So we are done
        break;
      }
    }
    trialHistory[[i]] <- sp23DesignObj$interimLookHistoryDF
    if (showProgress) setTxtProgressBar(pb, i)
  }
  if (showProgress) close(pb)
  trialHistory
}



###################################################
### code chunk number 1: Generate-Clinical-Trial-Data-Function
###################################################
generateClinicalTrialData <- function(nRec, nFUp, pi0, pi1,
                                      theta, lambda0, blockSize = 10) {
  ##
  ## nRec       is the number of patients recruited every year. Length(nRec) is the
  ##            number of years of recruitment
  ## nFUp       is the number of additional years of followup
  ## pi0         the probability of response under control arm
  ## pi1         the probability of response under treatment arm
  ## theta      the three dimensional parameter (alpha, beta, gamma) of the joint response/survival model
  ## lambda0    the baseline hazard rate
  ## blockSize  the size of the blocks for randomization of the treatment/control
  ##
  ## Returns a data frame with rows in order of patient entry
  ##

  N <- sum(nRec) # Total number recruited
  ##
  ## Check that blockSize is even and that blockSize divides N
  ##
  nBlocks <- N / blockSize
  if (trunc(blockSize / 2) * 2 != blockSize
      || trunc(nBlocks) * blockSize != N) {
    stop("Improper blocksize or number of subjects recruited")
  }

  ## Generate uniform entry times in each recruitment year
  entryTime <- sort(unlist(lapply(1:length(nRec), function(j) runif(nRec[j], j-1, j))))

  ## Generate treatment indicators (1 = treatment, 0 = control)
  treatmentIndicator <-
    unlist(lapply(1:nBlocks, function(x) sample(c(rep(0, blockSize / 2), rep(1, blockSize / 2)))))

  ## Generate Bernoulli responses with appropriate probabilities
  responseIndicator <- integer(N)
  responseIndicator[treatmentIndicator == 1] <- rbinom(sum(treatmentIndicator), 1, pi1)
  responseIndicator[treatmentIndicator == 0] <- rbinom(N-sum(treatmentIndicator), 1, pi0)

  ## Compute hazard rate per model
  rate = lambda0 *
    exp(theta$alpha * responseIndicator +
        theta$beta *  treatmentIndicator +
        theta$gamma * responseIndicator * treatmentIndicator)
  ## Generate time to event
  timeToEvent <- rexp(N, rate = rate)

  ## Construct data frame of result, name variables appropriately
  result <- data.frame(entryTime=entryTime, responseIndicator=responseIndicator,
                       treatmentIndicator=treatmentIndicator, timeToEvent=timeToEvent)
  rownames(result) <- paste("Pat", 1:N, sep=".")
  result
}


###################################################
### code chunk number 2: Generate-Interim-Data-Function
###################################################
generateInterimData <- function(clinicalTrialDF, interimTime, administrativeCensoringTime) {
  ##
  ## Given a clinical trial data frame (with variables
  ## entryTime, responseIndicator, treatmentIndicator, timeToEvent)
  ## and a calendar interimTime, generate a data frame that represents the interim
  ## data at the interimTime.
  ## administrativeCensoringTime is the time when the study ends
  ## Return data frame with timeToEvent calculated as if we at interimTime,
  ## and add a variable called eventTime, which is calendar time of event,
  ## that is, timeToEvent + entryTime
  ##
  ## Narrow to data of interest
  ##d <- clinicalTrialDF[clinicalTrialDF$entryTime <= interimTime, ]
  ##N <- nrow(d)
  d <- clinicalTrialDF
  ## perform a parallel minimum
  interimTimeToEvent <- pmin(d$timeToEvent,
                             administrativeCensoringTime,
                             pmax(0, interimTime - d$entryTime))  ## T_i(t)
  d$delta <- ifelse(interimTimeToEvent == d$timeToEvent, 1, 0) ## delta_i(t)
  d$timeToEvent <- interimTimeToEvent
  d$eventTime <- d$timeToEvent + d$entryTime ## calendar event time
  d[d$timeToEvent > 0, ] ##restrict to patients seen up to time t
}


###################################################
### code chunk number 11: Hessian
###################################################
hessian <- function(theta, pi, interimData) {
  m = nrow(interimData)
  x <- cbind(interimData$responseIndicator,
             interimData$treatmentIndicator,
             interimData$responseIndicator * interimData$treatmentIndicator)

  w1 <- theta[1] * x[, 1] + theta[2] * x[, 2] + theta[3] * x[, 3]
  expW1 <- exp(w1)
  w2 <- rev(cumsum(rev(expW1)))
  u <- matrix(0, nrow=m, ncol=3)
  for (jj in 1:3)
    for (ii in 1:m)
      u[ii, jj] <- sum(x[ii:m, jj] * expW1[ii:m]) / w2[ii]

  likdot <- t(interimData$delta) %*% (x - u)

  likddot <- matrix(0, nrow=3 ,ncol=3)
  u2 <- matrix(0, nrow=3, ncol=3)
  for (ii in 1:m) {
    u2[1:3 , 1:3] <- 0  ## zero out u2
    for (jj in ii:m) {
      u2 <- u2 + (x[jj, ] %*% t(x[jj, ])) * expW1[jj]
    }
    likddot <- likddot + interimData$delta[ii] * (u2 / w2[ii] - u[ii,] %*% t(u[ii,]))
  }
  likddot = -likddot

  pi0 <- pi[1]; pi1 <- pi[2];
  xi <- exp(theta)
  d0 <- (pi0*xi[1]+1-pi0)-(pi1*xi[1]*xi[3]+1-pi1)*xi[2]

  B <- (pi0 * xi[1] + 1-pi0) - (1 - pi1) * xi[2] - d0

  A1 <- pi0 * xi[1] / B - 1
  A2 <- -(1 - pi1) * xi[2] / B - 1
  A3 <- -1 / B
  D2 <- rbind(c(1, 0, 0), c(0, 1, 0), c(A1, A2, A3))

  D <- matrix(0, 3, 3)
  D[1, 1] <- -A1 * (A1 + 1)
  D[2, 2] <- -A2 * (A2 + 1)
  D[1, 2] <- D[2, 1] <- -(A1 + 1) * (A2 + 1)
  D[2, 3] <- D[3, 2] <- -A3 * (A2 + 1)
  D[1, 3] <- D[3, 1] <- -A3 * (A1 + 1)
  D[3, 3] <- -A3^2
  - (likdot[3] * D + t(D2) %*% likddot %*% D2)
}


###################################################
### code chunk number 8: Partial-Likelihood
###################################################
loglik2 <- function(theta, interimData) {
  w <- theta[1] * interimData$responseIndicator + theta[2] * interimData$treatmentIndicator +
    theta[3] * interimData$responseIndicator * interimData$treatmentIndicator ## alpha*Y + beta * Z + gamma * Y * Z
  w1 = rev(cumsum(rev(exp(w))))
  sum(interimData$delta * (w - log(w1)))
}


###################################################
### code chunk number 9: Partial-Likelihood
###################################################
loglik2.repar0 <- function(xi, interimData, pi0, pi1, eta.hyp = 0) {
  theta <- c(log(xi[1]), log(xi[2]),
  log((pi0 * xi[1] + 1 - pi0 - (1 - pi1) * xi[2] - eta.hyp) / (pi1 * xi[1] * xi[2])))
  w <- theta[1] * interimData$responseIndicator + theta[2] * interimData$treatmentIndicator +
    theta[3] * interimData$responseIndicator * interimData$treatmentIndicator ## alpha*Y + beta * Z + gamma * Y * Z
  w1 = rev(cumsum(rev(exp(w))))
  sum(interimData$delta * (w - log(w1)))
}


###################################################
### code chunk number 14: PerformInterimLookFunction
###################################################
performInterimLook <-
function(k, trueParameters, trialParameters, glrBoundary,
                               interimData, interimLookHistoryDF, argRejectH0R) {
#  cat("k = ", k, "\n")
  result <- c(rep(NA, ncol(interimLookHistoryDF)-4), FALSE, FALSE, FALSE, FALSE, NA)
  names(result) <- c(names(interimLookHistoryDF), "b.metas.Last")
  numLooks <- length(trialParameters$interimLookTime)
  pdiff.hyp <- trueParameters$pdiffHyp
  m <- nrow(interimData)
  interimNumberOfEvents = sum(interimData$delta)
  ## Compute some quantities
  responseSummary <- computeResponseSummary(interimData)
  piHat <- responseSummary[c("controlRespProp", "treatmentRespProp")]
  pi0Hat <- piHat[1] ; pi1Hat <- piHat[2]
  piHatH0 <- rep(responseSummary["pooledProp"], 2)
  w1 <- loglik1(piVec = piHat, responseSummary)
  ## Compute GLR for Response
  glrRespH0 <- w1 - loglik1(piVec = piHatH0, responseSummary)
  ## Record values
  result["m0"] <- responseSummary["m0"]; result["m1"] <- responseSummary["m1"]
  result["y0"] <- responseSummary["y0"]; result["y1"] <- responseSummary["y1"]
  result["pi0Hat"] <- pi0Hat; result["pi1Hat"] <- pi1Hat
  result["pi0HatH0"] <- piHatH0[1]; result["pi1HatH0"] <- piHatH0[2]
  result["glrRespH0"] <- glrRespH0

  rejectH0R <- argRejectH0R ## Make a current copy

  ## Skip if you have already rejected H_0^R
  if (!rejectH0R) { # This comes from previous interim look!
    rejectH0R <- ((2 * glrRespH0 >= glrBoundary[k, "RespEfficacy"]^2) && (pi0Hat < pi1Hat)) ## Response efficacy
    ##result["rejectH0R"] <- rejectH0R
    ## Now you have not rejected H_0^R, check for futility
    if (!rejectH0R) { ## This comes from current interim look
      initialPStart <- max(piHat[1] - responseSummary["m1"] * pdiff.hyp / responseSummary["m"], 0.01)
      test <- optim(par = initialPStart,
                    fn = loglik1GivenDelta,
                    respSummary = responseSummary, delta = pdiff.hyp,
                    lower = 0.001, upper = 0.999 - pdiff.hyp,
                    method="L-BFGS-B")
      piHatH1 <- c(test$par, test$par + pdiff.hyp)
      glrRespH1 <- w1 - loglik1(piVec = piHatH1, responseSummary)
      ## record values
      result["pi0HatH1"] <- piHatH1[1]; result["pi1HatH1"] <- piHatH1[2]; result["glrRespH1"] <- glrRespH1
      ##futility test equation 19 of paper
      acceptH0R <- ((2 * glrRespH1 >= glrBoundary[k, "RespFutility"]^2) && (pi1Hat < pi0Hat + pdiff.hyp))
      result["acceptH0R"] <- acceptH0R
    }
  }
  result["rejectH0R"] <- rejectH0R
  ## If, by the penultimate look H_0^R is not rejected, give up
  if ((k == (numLooks - 1)) && !rejectH0R) {
    acceptH0R <- TRUE
    result["acceptH0R"] <- acceptH0R
  } else {
    ## You might have rejected H_0^R in a previous look OR in the current look
    ## so you have test on the flag again!
    if (rejectH0R) {
      rejectH0S <- acceptH0S <- FALSE ## initially false
      if (k < numLooks) { # for all but final looks
        test <- optim(par=c(0,0,0), fn=loglik2,
                      interimData=interimData, control=list(fnscale=-1))
        ## Store the MLE theta.hat
        thetaHat <- test$par
        xiHat <- exp(thetaHat)
        hazard = computeDGivenXi(piHat, xiHat)
##        cat("thetaHat=", thetaHat, "\n")
        test <- constrOptim(theta=c(1.5, 0.5),
                            f=loglik2.repar0,
                            grad=NULL,
                            ui=rbind(c(1, 0), c(0, 1), c(pi0Hat, pi1Hat - 1)),
                            ci=c(0, 0, pi0Hat - 1 + 0),
                            control=list(fnscale=-1),
                            interimData=interimData, pi0=pi0Hat, pi1=pi1Hat)
        ## Now compute the c back (= w)
        w <- solveForCGivenABD(piHat, test$par[1], test$par[2], d=0)
        thetaHatH0 <- log(c(test$par, w))
        w1 = loglik2(theta = thetaHat, interimData)
        glrSurvH0 = w1 - loglik2(theta = thetaHatH0, interimData)
        if (hazard < trueParameters$etaHyp) {
          test <- constrOptim(theta=c(3.0, 0.2),
                              f=loglik2.repar0,
                              grad=NULL,
                              ui=rbind(c(1,0), c(0,1), c(pi0Hat, pi1Hat-1)),
                              ci=c(0, 0, trueParameters$etaHyp - (1 - pi0Hat) + 0),
                              control=list(fnscale=-1),
                              interimData=interimData, pi0=pi0Hat, pi1=pi1Hat,
                              eta.hyp=trueParameters$etaHyp)
          w = solveForCGivenABD(piHat, test$par[1], test$par[2], trueParameters$etaHyp)
          thetaHatH1=log(c(test$par,w))
          glrSurvH1 = w1 - loglik2(theta = thetaHatH1, interimData)
        } else {
          thetaHatH1 <- rep(NA, 3)
          glrSurvH1 = NA
        }
        ## Record quantities
        result["glrSurvH0"] <- glrSurvH0; result["glrSurvH1"] <- glrSurvH1
        result["alphaHat"] <- thetaHat[1]; result["betaHat"] <- thetaHat[2]; result["gammaHat"] <- thetaHat[3];
        result["hazard"] <- hazard
        result["alphaHatH0"] <- thetaHatH0[1]; result["betaHatH0"] <- thetaHatH0[2]; result["gammaHatH0"] <- thetaHatH0[3];
        #print("before")
        result["alphaHatH1"] <- thetaHatH1[1]; result["betaHatH1"] <- thetaHatH1[2]; result["gammaHatH1"] <- thetaHatH1[3];
        #print("after")

        ## Start test for efficacy
        ## Now have to take into account that we may have too few events
        ## THIS MEANS: ENSURE THAT NO ONE SPECIFIES TRIAL MIN NO OF EVENTS
        ## TO BE > any interim # of events
        if (interimNumberOfEvents >= trialParameters$minimumNumberOfEvents) {
          ## CHECK CHECK CHECK for proper args theta and pi
          v <- computeGammaSubT(thetaHat = thetaHatH0, pi = piHat, interimData)
          result["v"] <- v
          if (v > 0) {
            usedVIndices <- which(interimLookHistoryDF$usedV > 0)
            nUV <- length(usedVIndices)
            if (nUV > 0) {
              ## Only record if it is large enough compared to the previous recorded value
              if (v >= interimLookHistoryDF$usedV[usedVIndices[nUV]] * (1 + trialParameters$minimumIncreaseInV))
                result["usedV"] <- v
            } else {
              result["usedV"] <- v
            }
          }
          rejectH0S <- ((2 * glrSurvH0 >= glrBoundary[k, "SurvEfficacy"]^2) && (hazard > 0))
          result["rejectH0S"] <- rejectH0S
        }
        ## End test for efficacy
        ## Start test for futility
        acceptH0S <- ((2 * glrSurvH1 >= glrBoundary[k, "SurvFutility"]^2) && (hazard < trueParameters$etaHyp))
        result["acceptH0S"] <- acceptH0S
        ## End test for futility
      } else { ## final look
        test <- optim(par=c(0,0,0), fn=loglik2,
                      interimData=interimData, control=list(fnscale=-1))
        ## Store the MLE theta.hat
        thetaHat <- test$par
        xiHat <- exp(thetaHat)
        hazard = computeDGivenXi(piHat, xiHat)
##        cat("thetaHat=", thetaHat, "\n")
        test <- constrOptim(theta=c(1.5,0.5),
                            f=loglik2.repar0,
                            grad=NULL,
                            ui=rbind(c(1, 0), c(0, 1), c(pi0Hat, pi1Hat - 1)),
                            ci=c(0, 0, pi0Hat - 1 + 0),
                            control=list(fnscale=-1),
                            interimData=interimData, pi0=pi0Hat, pi1=pi1Hat)
        ## Now compute the c back (= w)
        w <- solveForCGivenABD(piHat, test$par[1], test$par[2], d=0)
        thetaHatH0 <- log(c(test$par, w))
        w1 = loglik2(theta = thetaHat, interimData)
        glrSurvH0 = w1 - loglik2(theta = thetaHatH0, interimData)
        ## Record quantities
        result["glrSurvH0"] <- glrSurvH0
        result["alphaHat"] <- thetaHat[1]; result["betaHat"] <- thetaHat[2]; result["gammaHat"] <- thetaHat[3];
        result["hazard"] <- hazard
        result["alphaHatH0"] <- thetaHatH0[1]; result["betaHatH0"] <- thetaHatH0[2]; result["gammaHatH0"] <- thetaHatH0[3];
        v <- computeGammaSubT(thetaHat = thetaHatH0, pi = piHat, interimData)
        result["v"] <- result["usedV"] <- v
        usedVIndices <- which(interimLookHistoryDF$usedV > 0)
        nUV <- length(usedVIndices)
##        cat("before nUV", nUV, "nonZeroVIndices", nonZeroVIndices, "\n")
##        cat("V So far", interimLookHistoryDF$v, "\n")
        usedV <- { if (nUV > 0) c(interimLookHistoryDF$usedV[usedVIndices], v) else v }
        ## NEED TO CALCULATE C=B.METAS[K]
        nUV <- length(usedV)
##        cat("after nUV", nUV, "usedV", usedV, "\n")
        if (nUV == 1) {
          b.metas.Last <- qnorm(1 - trialParameters$type1Error)
        } else if (nUV >= 2) {
          if (usedV[nUV] < usedV[nUV-1]){
            b.metas.Last <- Inf
          } else {
            test <- mHP.c(mu=numeric(nUV),
                          v=usedV,
                          b=glrBoundary[1, "SurvEfficacy"], alpha=trialParameters$type1Error,
                          eps=trialParameters$epsType1)
            b.metas.Last <- test$c
          }
        }
        result["b.metas.Last"] <- b.metas.Last
        ##        cat("c=", b.metas.Last,"\n")
        ##        cat("k=5 (final look): UsedV=",usedV,"\n")
        ## Check boundary crossing for c=b.metas[k]
        rejectH0S <- ((2 * glrSurvH0 >= b.metas.Last^2) && (hazard > 0))
        acceptH0S <- !rejectH0S
        result["rejectH0S"] <- rejectH0S
        result["acceptH0S"] <- acceptH0S
      }
    }
  }
  ## Return result vector
  result
}


###################################################
### code chunk number 16: Reset-Design-Function
###################################################
resetSP23Design <- function(sp23Design) {
  numLooks <- length(sp23Design$trialParameters$interimLookTime)
  naVector <- rep(NA, numLooks) # a qty we reuse
  logicalVector <- logical(numLooks) # another qty we reuse.

  interimLookHistoryDF <- data.frame(m0 = naVector,
                                     m1 = naVector,
                                     y0 = naVector,
                                     y1 = naVector,
                                     pi0Hat = naVector,
                                     pi1Hat = naVector,
                                     pi0HatH0 = naVector,
                                     pi1HatH0 = naVector,
                                     pi0HatH1 = naVector,
                                     pi1HatH1 = naVector,
                                     glrRespH0 = naVector,
                                     glrRespH1 = naVector,
                                     hazard = naVector,
                                     alphaHat = naVector,
                                     betaHat = naVector,
                                     gammaHat = naVector,
                                     alphaHatH0 = naVector,
                                     betaHatH0 = naVector,
                                     gammaHatH0 = naVector,
                                     alphaHatH1 = naVector,
                                     betaHatH1 = naVector,
                                     gammaHatH1 = naVector,
                                     glrSurvH0 = naVector,
                                     glrSurvH1 = naVector,
                                     v = naVector,
                                     usedV = naVector,
                                     rejectH0R =  logicalVector,
                                     acceptH0R =  logicalVector,
                                     rejectH0S =  logicalVector,
                                     acceptH0S =  logicalVector
                                     )
  list(trueParameters = sp23Design$trueParameters, trialParameters=sp23Design$trialParameters,
       glrBoundary=sp23Design$glrBoundary, interimLookHistoryDF=interimLookHistoryDF)
}


###################################################
### code chunk number 4: SolveForCGivenABD
###################################################
solveForCGivenABD <- function(piVec, a, b, d) {
  pi0 <- piVec[1]; pi1 <- piVec[2]
  ((pi0 * a + 1 - pi0) - (1 - pi1) * b - d) / (pi1 * a *b)
}


###################################################
### code chunk number 6: likelihood1
###################################################
loglik1 <- function(piVec, respSummary) {
  y0 <- respSummary["y0"]
  y1 <- respSummary["y1"]
  m0 <- respSummary["m0"]
  m1 <- respSummary["m1"]
  pi0 <- piVec[1]
  pi1 <- piVec[2]

  if ( ( (pi0 <= 0) && (y0 == 0) ) || ( (m0 == y0) && (pi0 >= 1) ) ) {
    term1PlusTerm2 <- 0
  } else {
    term1PlusTerm2 <- y0 * log(pi0) + (m0 - y0) * log(1 - pi0)
  }

  if ( ( (pi1 <= 0) && (y1 == 0) ) || ( (m1 == y1) && (pi1 >= 1) ) ) {
    term3PlusTerm4 <- 0
  } else {
    term3PlusTerm4 <- y1 * log(pi1) + (m1 - y1) * log(1 - pi1)
  }
  term1PlusTerm2 + term3PlusTerm4
}


###################################################
### code chunk number 7: maxlikelihood1
###################################################
loglik1GivenDelta <- function(p, respSummary, delta=0) {
  y0 <- respSummary["y0"]
  y1 <- respSummary["y1"]
  m0 <- respSummary["m0"]
  m1 <- respSummary["m1"]
  p1 <- p + delta

  if ( ( (p <= 0) && (y0 == 0) ) || ( (m0 == y0) && (p >= 1) ) ) {
    term1PlusTerm2 <- 0
  } else {
    term1PlusTerm2 <- -y0 * log(p) - (m0 - y0) * log(1 - p)
  }

  if ( ( (p1 <= 0) && (y1 == 0) ) || ( (m1 == y1) && (p1 >= 1) ) ) {
    term3PlusTerm4 <- 0
  } else {
    term3PlusTerm4 <- -y1 * log(p1) - (m1 - y1) * log(1 - p1)
  }
  term1PlusTerm2 + term3PlusTerm4
}


