CV = function(data, learner, params, fold=5, verbose=TRUE) {
  stopifnot(class(learner) == "CVST.learner" &&
            class(data) == "CVST.data" &&
            class(params) == "CVST.params")
  nParams = length(params)
  dimnames = list(as.character(1:fold), names(params))

  results = matrix(0, fold, nParams, dimnames=dimnames)
  size = getN(data) / fold

  for (ind in 1:nParams) {
    p = params[[ind]]
    for (f in 1:fold) {
      validationIndex = seq((f-1)*size + 1, f*size)
      curTrain = getSubset(data, -validationIndex)
      curTest = getSubset(data, validationIndex)
      # either mean squared error or mean classification error
      results[f, ind] = mean(.getResult(curTrain, curTest, learner, p))
    }
    if (verbose) {
      cat(names(params)[ind], "(", mean(results[, ind]), ")\n")
    }
  }
  winner = which.min(apply(results, 2, mean))
  if (length(winner) == 0) {
    return(NULL)
  }
  else {
    return(params[winner])
  }
}

# the function to perform fastcrossvalidation:
#
# train: training data CVST.data
#
# learner: the learner as CVST.learner
#
# params: list of parameters for the learner as CVST.params
#
# setup: setup of the CVST as CVST.setup
#
# test: either the test data for fixed test error setting or NULL, if
# the adjusted test error setting should be used
fastCV = function(train, learner, params, setup, test=NULL, verbose=TRUE) {
  stopifnot(class(learner) == "CVST.learner" && class(train) == "CVST.data" &&
            class(params) == "CVST.params" && class(setup) == "CVST.setup" &&
            (is.null(test) || class(test) == "CVST.data"))
  isClassificationTask = isClassification(train)
  regressionSimilarityViaOutliers = setup$regressionSimilarityViaOutliers
  earlyStopping = setup$earlyStoppingSignificance
  similarity = setup$similaritySignificance
  # use nested modeling, i.e. we start with the first minimalModel number of
  # data points and in each step subsequently add minimalModel data points to it
  nestModel = TRUE
  earlyStoppingWindow = setup$earlyStoppingWindow

  if (is.null(test)) {
    # we are in the adjusted test error setting, therefore we have to keep
    # an additional slice of the data for the last test
    minimalModel = getN(train) / (setup$steps + 1)
    n = getN(train) - minimalModel
  }
  else {
    minimalModel = getN(train) / setup$steps
    n = getN(train)
  }

  N = seq(minimalModel, n, by=minimalModel)
  st = getCVSTTest(setup$steps, setup$beta, setup$alpha)
  nParams = length(params)
  if (verbose) {
    cat("Total number of params:", nParams, "\n")
  }
  dimnames = list(names(params), as.character(N))
  traces = matrix(0, nParams, length(N), dimnames=dimnames)
  success = matrix(0, nParams, length(N), dimnames=dimnames)
  skipCalculation = rep(FALSE, nParams)
  isEarlyStopping = FALSE
  stoppedAt = length(N)
  activeConfigurations = matrix(FALSE, nParams, length(N), dimnames=dimnames)
  configurationsLeft = nParams
  
  for (ind in 1:length(N)) {
    n = N[ind]
    if (!isClassificationTask && regressionSimilarityViaOutliers) {
      err = .calculateErrors(train, test, n, learner, params, skipCalculation, squared=FALSE)
      success[, ind] = apply(err^2, 1, mean)
    }
    else {
      err = .calculateErrors(train, test, n, learner, params, skipCalculation)
      success[, ind] = apply(err, 1, mean)
    }
    
    success[, ind] = apply(err, 1, mean)
    indByError = sort.list(success[, ind], decreasing=FALSE, na.last=TRUE)
    traces[indByError[1], ind] = 1
    sortedErrors = t(err[indByError, ])
    if (!isClassificationTask && regressionSimilarityViaOutliers) {
      s = apply(sortedErrors, 2, sd)
      sortedErrors = t(abs(t(sortedErrors)) > s * qnorm(1 - (similarity / 2)))
    }
    adjustedSignificance = similarity / (configurationsLeft - 1)
    for (k in 2:length(indByError)) {
      if (is.na(success[indByError[k], ind])) {
        # we either have an unsufficient model, which gives us NA as result...
        # ... or reached the skipCalculation, so we can stop our procedure
        break
      }
      if (isClassificationTask) {
        pvalue = cochranq.test(sortedErrors[, 1:k])$p.value
      }
      else {
        if (regressionSimilarityViaOutliers) {
          pvalue = cochranq.test(sortedErrors[, 1:k])$p.value
        }
        else {
          pvalue = friedman.test(sortedErrors[, 1:k])$p.value
        }
      }
      if (!is.nan(pvalue) && pvalue <= adjustedSignificance) {
        break
      }
      traces[indByError[k], ind] = 1
    }
    if (verbose) {
      cat("(sim:", sum(traces[, ind]), "alpha:", similarity, "left:", configurationsLeft, ")")
    }
    # do the testing here...
    # check for loosers
    if (ind > 1) {
      testResults = apply(traces[, 1:ind], 1, testSequence, st=st)
      # check for loosers
      skipCalculation = (testResults == -1)
      if (verbose) {
        cat("Skipped configurations:", sum(skipCalculation), " ")
      }
    }
    configurationsLeft = nParams - sum(skipCalculation)
    activeConfigurations[, ind] = !skipCalculation
    # check for early stopping
    if (earlyStoppingWindow >= 2 && ind > earlyStoppingWindow && earlyStopping < 1.0) {
      # check, whether all remaining parameters perform similar
      if (sum(!skipCalculation) > 1)
        pvalue = cochranq.test(t(traces[!skipCalculation, (ind-earlyStoppingWindow+1):ind]))$p.value
      else {
        pvalue = 1.0
      }
      if (!is.nan(pvalue) && pvalue > earlyStopping) {
        if (verbose) {
          cat("EARLY STOPPING!")
        }
        isEarlyStopping = TRUE
        stoppedAt = ind
        break
      }
      # just go on, if they are signifcantly dissimilar!
    }
  }
  if (verbose) {
    cat("\n")
  }
  theWinners = !skipCalculation
  ret = list(traces=traces, success=success)
  ret$numberOfPotentialWinners = sum(theWinners)
  ret$isEarlyStopping = isEarlyStopping
  ret$stoppedAt = stoppedAt
  ret$activeConfigurations = activeConfigurations
  ret$earlyStoppingWindow = earlyStoppingWindow
  winningConfiguration = .getOptimalSolution(ret)

  ret$param = params[winningConfiguration]
  ret$winningConfiguration = winningConfiguration

  return(params[winningConfiguration])
}

# returns a (# configuration) X (# testsamples) matrix containing 0/1 or squared error at
# position i, j if the model learned on N data points of traindata
# with configuration i labeled point j of the testdata
# correctly or not. skipCalculation controls, which confguration should be
# skipped. A NA in the returned matrix corresponds to skipped configuration.
.calculateErrors = function(traindata, testdata, N, learner, params, skipCalculation, squared=TRUE) {
  nestModel = TRUE
  nPars = length(params)

  if (nestModel) {
    sampleIndex = 1:N
  }
  else {
    sampleIndex = sample.int(getN(traindata), N)
  }
  # if no test data is available, we have the adjusted test error settings,
  # i.e. we use the rest of the train data, which is not used for model building
  # to determine the test error
  if (is.null(testdata)) {
    testdata = getSubset(traindata, -sampleIndex)
  }
  # initialize results
  results = matrix(NA, nPars, getN(testdata))
  # calculate results
  curTrain = getSubset(traindata, sampleIndex)
  for (ind in 1:nPars) {
    param = params[[ind]]
    if (!is.null(skipCalculation) && skipCalculation[ind]) {
      next
    }
    results[ind, ] = as.vector(.getResult(curTrain, testdata, learner, param, squared=squared))
  }
  return(results)
}

.getOptimalSolution = function(paramRace) {
  remainingConfs = paramRace$activeConfigurations[, paramRace$stoppedAt]
  if (sum(remainingConfs) == 1) {
    return(remainingConfs)
  }
  # pick the race, which has the smallest mean rank inside
  # the earlyStoppingWindow:
  lastSuccess = paramRace$success[remainingConfs, (paramRace$stoppedAt - paramRace$earlyStoppingWindow + 1):paramRace$stoppedAt]
  meanRank = apply(apply(lastSuccess, 2, rank), 1, mean)
  # breaks ties at random
  overallWinner = which(remainingConfs)[.which.is.min(meanRank)]

  ret = rep(FALSE, nrow(paramRace$traces))
  names(ret) = rownames(paramRace$traces)
  ret[overallWinner] = TRUE
  return(ret)
}

.which.is.min = function (x) {
  y = seq_along(x)[x == min(x)]
  if (length(y) > 1) {
    y = sample(y, 1)
  }
  return(y)
}

getCVSTTest = function(steps, beta=.1, alpha=.01) {
  pi1 = .5 * ((1 - beta) / alpha)^(1/steps)
  sst = constructSequentialTest(.5, pi1, beta, alpha)
  sst$steps = steps
  return(sst)
}
