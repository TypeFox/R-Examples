# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# These helper functions check parameter integrity

.checkIntParam = function(param, paramName, positive)
{
  if ((!missing(positive) && param < (if (positive) 1 else 0)) ||
      (param != floor(param)) ||
      (length(param) > 1))
  {
    if (missing(positive))  
    {
      paste("\n*", paramName, "must contain a single integer instead of", param)
    }
    else if (positive)
    {
      paste("\n*", paramName, "must contain a single positive integer instead of", param)
    }
    else 
    {
      paste("\n*", paramName, "must contain a single non-negative integer instead of", param)
    }
  }
}


.checkLogicalParam = function(param, paramName)
{
  if (length(param) > 1 ||
      !is.logical(param))
  {
    paste("\n* ", paramName, " must contain a single logical value (i.e., TRUE or FALSE)", sep="")
  }
}


# This helper function returns the probabilities of each element of eventList

.getEventListProbs = function(ndicePerRoll, nsidesPerDie, eventList)
{
  probs = getSumProbs(ndicePerRoll, nsidesPerDie)$probabilities

  # On the assumption that eventList has length nrolls (which is safe since this is a 
  # private helper function), we calculate the probability of getting an acceptable 
  # outcome (a "success") on each of the rolls by iterating through the vector of 
  # successes for that roll and adding the corresponding probability to our tally

  eventListProbs = c()
  for (i in 1:length(eventList)) 
  {
    successesForThisRoll = sort(eventList[[i]])
    successProbForThisRoll = 0
    for (j in 1:length(successesForThisRoll))
    {
      successProbForThisRoll = successProbForThisRoll + probs[(successesForThisRoll[j] - (ndicePerRoll - 1)),2]
    }
    eventListProbs[i] = successProbForThisRoll
  }
  eventListProbs 
}


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# NOTE: the parameters nrolls, ndicePerRoll, nsidesPerDie, and nkept in the function
# signatures below depart slightly from our usual coding conventions (i.e., are not 
# numRolls, numDicePerRoll, numSidesPerDie, and numDiceKept) so that they are easily 
# abbreviated as "nr", "nd", "ns", and "nk", respectively, in function calls

getEventProb = function(nrolls,
                        ndicePerRoll,
                        nsidesPerDie,
                        eventList,
                        orderMatters=FALSE)
{
  errorVector = character()
  errorVector = append(errorVector, .checkIntParam(nrolls, "nrolls", positive=TRUE))
  errorVector = append(errorVector, .checkIntParam(ndicePerRoll, "ndicePerRoll", positive=TRUE))
  errorVector = append(errorVector, .checkIntParam(nsidesPerDie, "nsidesPerDie", positive=TRUE))

  if (length(eventList) > nrolls)
  {
    errorVector = append(errorVector, "\n* The length of eventList must not be greater than nrolls")
  }
  if (orderMatters & length(eventList) != nrolls)
  {
    errorVector = append(errorVector, "\n* If orderMatters is passed as TRUE, the length of eventList must equal\n  nrolls (i.e., there must be an element of eventList for each roll)")
  }
  if (!all(sapply(eventList, is.numeric)))
  {
    errorVector = append(errorVector, "\n* All elements of eventList must be numeric vectors")
  }
  if (!all(as.logical(sapply(sapply(eventList, function(x) x == floor(x)), min))))
  {
    errorVector = append(errorVector, "\n* All numbers in each element of eventList must be positive integers")
  }
  if (min(sapply(eventList, min)) < ndicePerRoll ||
      max(sapply(eventList, max)) > (ndicePerRoll * nsidesPerDie))
  {
    errorVector = append(errorVector, "\n* All numbers in each element of eventList must be between ndicePerRoll\n  and (ndicePerRoll * nsidesPerDie)")
  }
  errorVector = append(errorVector, .checkLogicalParam(orderMatters, "orderMatters"))

  if (length(errorVector) > 0)
  {
    stop(errorVector)
  }

  eventList = lapply(eventList, unique)

  # If eventList doesn't have an element for each roll, we add elements until it does;
  # after this point, each element of eventList will constrain one roll (but some of 
  # those constraints may be simply {min:max} for that roll--i.e., trivial constraints)

  if (length(eventList) < nrolls)
  {
    eventList = lapply(c(eventList, rep(0, nrolls - length(eventList))), function(x){x = (if (max(x) == 0) ndicePerRoll:(ndicePerRoll*nsidesPerDie) else x)})    
  }

  if (orderMatters)
  {
    outcomeProb = prod(.getEventListProbs(ndicePerRoll, nsidesPerDie, eventList))
  }
  else # i.e., if (!orderMatters)
  {
    # We only calculate probabilities if each element of eventList is a length-1 vector
    # (i.e., a single number), e.g., {2, 3, 2}; if any element is longer than that, e.g.,
    # {2, {3, 4}, 2}, we call ourselves recursively on each list we can construct of only
    # length-1 vectors (e.g., in the example above we'd call ourselves on {2, 3, 2} and 
    # {2, 4, 2}); then we sum the resulting probabilities (which, since orderMatters is 
    # FALSE, account for all permutations of each of {2, 3, 2} and {2, 4, 2}) to arrive 
    # at our probability for the original list of {2, {3, 4}, 2}

    listElemLengths = sapply(eventList, length)
    maxListElemLength = max(listElemLengths)
    if (maxListElemLength > 1)
    {
      # Here we populate combMatrix with the elements of eventList to produce a 
      # matrix each row of which is a selection of one element from each element of
      # eventList; e.g., given the eventList {{1, 2}, {1, 2, 4}, 2}, we'd produce
      # a 6 x 3 matrix with rows {1, 1, 2}, {1, 2, 2}, {1, 4, 2}, {2, 1, 2}, {2, 2, 2},
      # and {2, 4, 2}

      combMatrix = matrix(nrow = prod(listElemLengths), ncol = nrolls)
      if (nrolls > 1)
      {
        for (i in 1:(nrolls-1))
        {
          combMatrix[,i] = rep(eventList[[i]], each = prod(listElemLengths[(i+1):nrolls]))
        }
      }
      combMatrix[,nrolls] = rep(eventList[[nrolls]])

      # Next we eliminate all rows that are permutations of other rows (otherwise we
      # would over-count in the calculations that follow)

      if (nrolls > 1)
      {
        combMatrix = unique(t(apply(combMatrix,1,sort)))
      }
      else
      {
        combMatrix = unique(combMatrix)
      }

      # Now we make a recursive call for each row of combMatrix and sum the resulting
      # probabilities to arrive at our probability for the original eventList

      sumOfProbs = sum(apply(combMatrix,
                             1,
                             function(x) getEventProb(nrolls,
                                                      ndicePerRoll,
                                                      nsidesPerDie,
                                                      as.list(x),
                                                      orderMatters)))
      outcomeProb = sumOfProbs
    }
    else
    {

      # If each element of eventList is a length-1 vector, we can convert eventList
      # itself to a vector; then we calculate the probability of getting the specific
      # set of outcomes specified by eventList in any order (reflecting the fact that 
      # orderMatters was passed in as FALSE)
     
      eventListAsVector = sapply(eventList, max)
      eventListProb = prod(.getEventListProbs(ndicePerRoll, nsidesPerDie, eventListAsVector))
      outcomeProb = eventListProb * factorial(nrolls) / prod(factorial(table(eventListAsVector)))
    }
  }

  outcomeProb
}


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


getSumProbs = function(ndicePerRoll,
                       nsidesPerDie,
                       nkept = ndicePerRoll,
                       dropLowest = TRUE,
                       sumModifier = 0,
                       perDieModifier = 0,
                       perDieMinOfOne = TRUE)
{

  # We begin with preliminary error-checking

  errorVector = vector(mode = "character", length = 0)
  errorVector = append(errorVector, .checkIntParam(ndicePerRoll, "ndicePerRoll", positive=TRUE))
  errorVector = append(errorVector, .checkIntParam(nsidesPerDie, "nsidesPerDie", positive=TRUE))
  errorVector = append(errorVector, .checkIntParam(nkept, "nkept", positive=TRUE))
  if (nkept > ndicePerRoll)
  {
    errorVector = append(errorVector, "\n* nkept must not be greater than ndicePerRoll")
  }
  errorVector = append(errorVector, .checkIntParam(sumModifier, "sumModifier"))
  errorVector = append(errorVector, .checkIntParam(perDieModifier, "perDieModifier"))
  errorVector = append(errorVector, .checkLogicalParam(perDieMinOfOne, "perDieMinOfOne"))

  if (length(errorVector) > 0)
  {
    stop(errorVector)
  }

  numOutcomes = nsidesPerDie^ndicePerRoll
  numDiceToDrop = ndicePerRoll - nkept
  currNumArrangements = 0
  
  sumModifier = sumModifier + (perDieModifier * nkept)
  
  currentSum = 0

  vectorOfSums = as.integer((nkept + sumModifier) :
                                        ((nsidesPerDie * nkept) + sumModifier))

  numPossibleSums = length(vectorOfSums)


  # sumTallyMatrix is used to track the number of times we see every possible outcome sum,
  # which we will use to produce the probabilities of every sum (e.g., for the 3d6 case we 
  # see 10 as a sum 27 times, so the probability of a sum of 10 is 27/216 = .125, while for
  # the 5d6 drop 2 case we see 13 as a sum 1055 times, so the probability of a sum of 13
  # is 1055/7776 = .1356739).

  sumTallyMatrix = matrix(data = c(vectorOfSums,
                                   as.integer(rep.int(0, numPossibleSums)),
                                   as.integer(rep.int(0, numPossibleSums))),
                          nrow = numPossibleSums,
                          ncol = 3,
                          dimnames = list(NULL,
                                          c("Sum",
                                            "Probability",
                                            "Ways to Roll")))

  # boundaryVal is the most extreme die-roll value that will be kept (i.e., the die-roll "boundary"
  # value: e.g., for 5d6 drop lowest two, if our sorted die rolls are {3 4 4 5 6}, boundaryVal is 4).
  # We'll call all dice with this value the "b" dice (because they're on the [b]oundary).

  for (boundaryVal in 1 : nsidesPerDie)
  {

    # numOs is the number of dice whose values are outside of boundaryVal (e.g., for 5d6 drop lowest 
    # two, if we roll {3 4 4 5 6}, boundaryVal is 4, so numOs is 1).  We'll call these dice the "o"
    # dice (because they're [o]utside our boundary).  
    # NOTE: We have an embedded if clause in our for-loop declaration because if we're dropping lowest 
    # and boundaryVal is 1 or we're dropping highest and boundaryVal is nsidesPerDie, there cannot be 
    # any dice whose values are outside of boundaryVal, and hence numOs can only be 0.  The following 
    # loop syntax might look suspicious, but we *do* want to iterate once in these two cases, as well
    # as in the case where numDiceToDrop is 0 (and in all three such cases, numOs will be 0).

    for (numOs in 0 : (if ((dropLowest && boundaryVal == 1) ||
                           (!dropLowest && boundaryVal == nsidesPerDie)) 0 else numDiceToDrop))
    {

      # numBsKept is the number of b's that will be kept (e.g., for 5d6 drop lowest two, if we roll
      # {3 4 4 5 6}, numBsKept is 1, because one of the two 4's will be kept).
      # Now, since we're discarding numDiceToDrop dice (including the numOs o's), we'll discard 
      # (numDiceToDrop - numOs) of the b's and keep numBsKept of them, and thus the total number of
      # b's is (numBsKept + numDiceToDrop - numOs).  NOTE: Hence, the number of dice whose values
      # exceed boundaryVal is (nkept - numBsKept).  We will call these higher dice the "i" dice (since
      # they're "inside" our boundary).

      for (numBsKept in 1 : nkept)
      {

        # By this part of the function, we've specified a class of outcomes identified by their 
        # (boundaryVal, numOs, numBsKept) values--i.e., every outcome in this class has the
        # following properties:
        # 1). the die-roll boundary value boundaryVal; 
        # 2). numOs "o" dice, whose values are outside boundaryVal and will be dropped; and
        # 3). numBsKept "b" dice that will be kept.  Furthermore, each such outcome has
        # 4). (numBsKept + numDiceToDrop - numOs) "b" dice in total, and 
        # 5). (nkept - numBsKept) "i" dice, whose values are inside boundaryVal.

        numBs = (numBsKept + numDiceToDrop - numOs)
        numIs = (nkept - numBsKept)

paste("\n\nIn this class, boundaryVal is ", boundaryVal, ", numOs is ", numOs,", numBsKept is", numBsKept, ", numBs is ", numBs, ", and numIs is ", numIs, "\n", sep="")

        # Now, we're interested in sums for the various outcomes in this class, and these sums
        # don't depend upon the order in which the various values appear in our sequence of
        # rolls; i.e., multiple outcomes in this class will have the same values but have the o's,
        # the b's, and the i's appear at different places in the sequence of rolls.  To account
        # for this, we need to multiply each distinct outcome by the number of other outcomes
        # that are identical to it except for the order in which the o's, b's, and i's occur
        # (NOTE: the orders *within* these groups are accounted for below: order within the o's 
        # is accounted for immediately below, and we account for the order within the i's in the
        # section of the code in which we enumerate the i's).  For now, we will define a term by
        # which to multiply each sum we find; this term is a result of the multinomial theorem's 
        # combinatoric interpretation as the number of ways to put n distinct objects (in this 
        # case, our die rolls) into 3 bins of size numOs, (numBsKept + numDiceToDrop - numOs), 
        # and (nkept - numBsKept), corresponding to the number of o's, b's, and i's in the class:

        numArrangementsOfDice = (factorial(ndicePerRoll) / 
                                 (factorial(numOs) * factorial(numBs) * factorial(numIs)))

        # [NOTE: The formula above could overflow if ndicePerRoll gets large, in which case we may
        #  consider using lfactorial()--but I think the function would keel over before then anyway]

        # Because we support dropping lowest or highest values, we define convenient variables
        # to allow us to operate over appropriate ranges for the rest of this function

		innerBoundary = if (dropLowest) nsidesPerDie else 1
		outerBoundary = if (dropLowest) 1 else nsidesPerDie
        rangeOfOVals = abs(boundaryVal - outerBoundary)
        rangeOfIVals = abs(boundaryVal - innerBoundary)
        possibleIValsVec = if (dropLowest) ((boundaryVal+1) : nsidesPerDie) else (1 : (boundaryVal-1))

        # Next: The value of boundaryVal is fixed for this loop, but there are many "o" dice values
        # that outcomes in this class might have; because we don't care about these values, we need
        # to increase our multiplicity term to account for the outcomes that will be identical 
        # to this iteration's distinct outcome but for the values (and order) of the o's:

        numArrangementsOfDice = numArrangementsOfDice * rangeOfOVals^numOs

        # Now that we've accounted for sorting the values into three bins and for all the possible
        # "o" values that are immaterial to our calculations, we can treat our outcome class as
        # sorted into groups and can focus our attention on the numBsKept b's we keep and the 
        # i's.  The numBsKept b's will contribute a known amount to our sum for all outcomes in 
        # this class (viz., numBsKept * boundaryVal); but the i's will contribute various amounts, 
        # depending on their values.  So now we turn to determining the possible distinct outcomes 
        # for this class by enumerating the possible values for the i's.  We will work as follows:
        # rangeOfIVals is the distance between the smallest and largest possible "i" values for this
        # class of outcomes, and we use it to determine the number of distinct outcomes for this
        # class, which is given by rangeOfIVals^numIs.  We create an outcomeMatrix with as many
        # rows as there are distinct outcomes for this class and nkept columns; each element
        # in a row corresponds to a die-roll value, and the sum of the row elements is the
        # sum for that distinct outcome.  We populate outcomeMatrix with a row for every possible 
        # value for the i's in this class (and hence all distinct outcomes in the class).  We then
        # calculate the number of permutations of each distinct outcome (e.g., in the 3d6 case,
        # the outcome {1, 1, 2} has three permutations) and use this information to calculate the 
        # probability of every possible outcome in this class.

        if (numBsKept == nkept)
        {
          currentSum = (numBsKept * boundaryVal) + sumModifier
          # We adjust row index by (nkept - 1) so that, e.g., a 3d6 tally starts in row 1, not 3
          sumTallyMatrix[currentSum - sumModifier - (nkept - 1), 2] = sumTallyMatrix[currentSum - sumModifier - (nkept - 1), 2] + numArrangementsOfDice
        }
        else
        {
          outcomeMatrix = matrix(nrow = choose((rangeOfIVals + numIs - 1), numIs), ncol = nkept)

          if (dim(outcomeMatrix)[1] > 0)
          {

            outcomeMatrix[,1 : numBsKept] = boundaryVal

            hCombs = combinations(n = rangeOfIVals,
                                  r = numIs,
                                  v = possibleIValsVec,
                                  repeats.allowed = TRUE)
            hPermCounts = apply(hCombs, 1, function(x) factorial(numIs)/prod(factorial(table(x))))

            outcomeMatrix[,(numBsKept+1) : nkept] = hCombs

            for (rowNum in 1 : nrow(outcomeMatrix))
            {
              currentSum = sum(outcomeMatrix[rowNum,]) + sumModifier
              currNumArrangements = numArrangementsOfDice * hPermCounts[rowNum]
              sumTallyMatrix[currentSum - sumModifier - (nkept - 1), 2] = sumTallyMatrix[currentSum - sumModifier - (nkept - 1), 2] + currNumArrangements
            }
          }
        }
      }
    }
  }

  if (perDieMinOfOne)
  {
    if (sumTallyMatrix[numPossibleSums,1] <= nkept)
    {
      sumTallyMatrix = matrix(data = c(nkept, numOutcomes, numOutcomes),
                              nrow = 1,
                              ncol = 3,
                              dimnames = list(NULL,
                                              c("  Sum  ","  Probability  ", "  Ways to Roll  ")))
    }
    else
    {
      extraWaysToRollMin = sum(sumTallyMatrix[sumTallyMatrix[,1] < nkept,2])
      sumTallyMatrix = sumTallyMatrix[sumTallyMatrix[,1] >= nkept,]
      sumTallyMatrix[1,2] = sumTallyMatrix[1,2] + extraWaysToRollMin
    }
  }

  sumTallyMatrix[,3] = sumTallyMatrix[,2]
  sumTallyMatrix[,2] = sumTallyMatrix[,2] / numOutcomes

  overallAverageSum = sum(sumTallyMatrix[,1] * sumTallyMatrix[,3] / numOutcomes)

  
  list(probabilities = sumTallyMatrix, average = overallAverageSum)
}
