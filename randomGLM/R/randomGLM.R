# Peter Langfelder's additions and changes

# Main change: unifying all levels of interactions in a single function. 

# Introducing the concept of an interaction matrix that specifies how the input
# features get multiplied into the predictors used in the underlying models. 
# The matrix has maxInteractionOrder rows and as many columns as needed. Each
# column corresponds to one term in for the (g)lm models. The entries in a
# column give indices of the input features that are multiplied in the term; 0
# index means "none". Thus, if one has 3 input features x1, x2, x3, an interaction
# matrix with up to 3-way interactions may look like (first line is column name)

# c1 c2 c3 c4 c5 c6 c7
#  1  2  3  1  1  1  2 
#  0  0  0  1  2  3  2
#  0  0  0  0  2  0  0

# This means that term 1 (column 1) is x1, term 2 is x2, term 3 is x3, term 4 is
# x1*x1, term5 is x1*x2*x2, term 6 in x1*x3, term 7 is x2*x2. No particular order
# is assumed in any of the rows or columns. The function that generates this
# matrix, named .interactionMatrix, 
# also assigns unique names to each column which is helpful when identifying
# individual terms and extracting their content.
#
# The advantage of this representation is that it makes it very convenient to
# actually generate the interaction terms from given numeric data. The generation
# of the interactions is implemented in function .generateInteractions.
#
# Furthermore, given an interaction matrix, it is easy to count the occurences of
# each input feature. This is implemented in function .countsInInteractionMatrix
# which counts the number of times each input feature appears in each level of
# interactions in a given interaction matrix.


# . Removed scaling of predictors

## change to 3.0
# 1. function name change to randomGLM
# 2. parameter change: corFncForCandidateCovariates, corOptionsForCandidateCovariates
# 3. change all "gene" to "feature"
# 4. add interceptOfForwardRegression

# For now, include required packages

#library(gtools)
#library(MASS)

#=====================================================
#
# Helper functions
#
#=====================================================

.spaste = function(...)
{
  paste(..., sep = "");
}

# prettier print than the standard "print" function.

.cat.nl = function(...)
{
  cat(.spaste(..., "\n"));
}

.prependZeros = function(x, len = max(nchar(x)))
{
  lengths = nchar(x);
  if (len < max(lengths)) stop("Some entries of 'x' are too long.");
  out = as.character(x);
  n = length(x);
  for (i in 1:n) if (lengths[i] < len)
    out[i] = .spaste( paste(rep("0", len-lengths[i]), collapse = ""),
                     x[i]);

  out;
}

  
#=====================================================
#
# Generator of interactions
#
#=====================================================

# My own function for generating combinations since the one in gtools fails with large n's.

# Combinations with repeats

.combinations = function(n, order)
{
  if (order==1) return(matrix(c(1:n), 1, n));

  # If order is not 1, calculate result using recursion
  nOut = choose(n+order-1, order);
  out = matrix(0, order, nOut);
  sub = .combinations(n, order-1)
  index = 1;
  for (i in 1:n)
  {
    n1 = ncol(sub);
    out[, index:(index + n1 - 1)] = rbind( rep(i, ncol(sub)), sub);
    index = index + n1;
    sub = sub[, colSums(sub==i) == 0, drop = FALSE];
  }
  out;
}
 

# Generates an interaction matrix of all interactions up to specified maxOrder of
# variables indexed 1:n. Optionally also sets column names with some flexibility.
# This uses function cobinations from gtools

.interactionMatrix = function(n, maxOrder, 
                 setColNames = TRUE,
                 originalNames = c(1:n),
                 featureSeparator = ".")
{
  out = NULL;
  for (o in 1:maxOrder)
  {
    combs = .combinations(n, o);
    nameMatrix = array(originalNames[combs], dim = dim(combs));
    colnames = apply(nameMatrix, 2, paste, collapse = featureSeparator);
    if (o < maxOrder)
       combs = rbind( combs, matrix(0, maxOrder-o, ncol(combs)));
    if (setColNames) colnames(combs) = colnames;
    out = cbind(out, combs);
  }
  rownames(out) = .spaste("Feature.", c(1:maxOrder));
  out
}

# utility function used by .generateInteractions below.
.generateInteractions.1row = function(x, interactionMatrix)
{
  mat = x[interactionMatrix];
  dim(mat) = dim(interactionMatrix);
  apply(mat, 2, prod);
}

# Generates interactions from a given matrix of numeric data for "original" or
# "input" features and an interaction matrix. 
# the trick to generate everything in one step is to add a dummy "feature" equal to
# 1 in all rows to the actual features and use it as "feature with index 0". 
# The function optionally transfers the column names from the interaction matrix to
# the result, which makes the columns names of the result unique and somewhat
# descriptive.

.generateInteractions = function(x, maxOrder, interactionMatrix = NULL, x1 = NULL, 
                                 setColNames = FALSE,
                                 originalNames= c(1:ncol(x)))
{
  n = ncol(x);
  if (is.null(interactionMatrix)) 
    interactionMatrix = .interactionMatrix(n, maxOrder, setColNames = setColNames, 
                                           originalNames = originalNames);

  if (nrow(interactionMatrix)!=maxOrder) stop("Internal error: nrow(interactionMatrix)!=maxOrder.");
  if (maxOrder==1)
  { 
    # If maxOrder is 1, the result is trivial. No need to spend time generating a trivial product.
    x.int = x[, as.vector(interactionMatrix), drop = FALSE];
  } else {
    # The non-trivial case
    if (is.null(x1)) x1 = cbind(x, rep(1, nrow(x)));
    interactionMatrix[ interactionMatrix==0 ] = n+1;

    # If the number of variables (columns in interactionMatrix) is 1, the result needs no transposing;
    # otherwise it needs to be transposed.
    if (ncol(interactionMatrix)==1)
    {
      x.int = as.matrix(apply(x1, 1, .generateInteractions.1row, interactionMatrix))
    } else 
      x.int = t(apply(x1, 1, .generateInteractions.1row, interactionMatrix));
  }
  if (setColNames) colnames(x.int) = colnames(interactionMatrix);
  x.int
}

# This function counts the number of times each input feature appears in each level
# of interactions in a given interaction matrix.

.countsInInteractionMatrix = function(im, nFeatures)
{
  maxLevel = nrow(im);
  counts = matrix(0, maxLevel, nFeatures);
  level = maxLevel - colSums(im==0);
  for (l in 1:maxLevel)
  {
    mat1 = im[ 1:l, level==l, drop = FALSE];
    mat1.unique = sapply(as.data.frame(mat1), unique);
    counts1 = table(unlist(mat1.unique));
    where = as.numeric(names(counts1));
    counts[ l, where] = as.numeric(counts1);
  }
  rownames(counts) = .spaste("Level.", c(1:maxLevel));
  counts;
}


# Translate ordinal data using a dictionary

.translate = function(data, dictionary)
{
  translated = dictionary[ match(data, dictionary[, 1]), 2];
  attributes(translated) = attributes(data);
  translated;
}

#=======================================================================================================
#
# binarizeCategoricalVar
#
#=======================================================================================================
# Assumes x is a vector but can easily be modified to also work with matrices.
.binarizeCategoricalVar = function(x, minCount = 3, val1 = 0, val2 = 1, nameSep = ".vs.", namePrefix = "",
                                  nameForAll = "all", 
                                  ignore = NULL, includePairwise = TRUE,
                                  includeLevelVsAll = FALSE, levelOrder = NULL)
{
  tab = table(x);
  levels0 = names(tab);
  tab = tab[ tab >= minCount & !(levels0 %in% ignore) ];
  levels = names(tab);
  if (!is.null(levelOrder))
  {
    order = match(levelOrder, levels);
    order = order[is.finite(order)];
    levels0 = levels[order];
    levels1 = levels[ !levels %in% levels0];
    levels = c(levels0, levels1);
  }
  nSamples = length(x);
  nLevels = length(levels)
  nBinaryVars = includePairwise * nLevels * (nLevels - 1)/2 + includeLevelVsAll * nLevels;
  if (nBinaryVars==0) return(NULL);
  out = matrix(NA, nSamples, nBinaryVars)
  levelTable = matrix("", 2, nBinaryVars);
  ind = 1;
  names = rep("", nBinaryVars);
  if (includePairwise)
  {
    for (v1 in 1:(nLevels-1)) for (v2 in (v1+1):nLevels)
    {
       out[ x==levels[v1], ind] = val1;
       out[ x==levels[v2], ind] = val2;
       names[ind] = .spaste(namePrefix, levels[v1], nameSep, levels[v2]);
       levelTable[, ind] = levels[ c(v1, v2)];
       ind = ind + 1;
    }
  }
  if (includeLevelVsAll)
    for (v1 in 1:nLevels)
    {
      out[, ind] = c(val1, val2) [ as.numeric(x==levels[v1])+1 ];
      names[ind] = .spaste(namePrefix, nameForAll, nameSep, levels[v1]);
      levelTable[, ind] = c(nameForAll, levels[v1]);
      ind = ind+1;
    }
  colnames(out) = names;
  colnames(levelTable) = names;
  rownames(levelTable) = .spaste("Value.", c(val1, val2));
  attr(out, "includedLevels") = levelTable;
  out;
}



#=====================================================
#
# forwardSelection
#
#=====================================================

.forwardSelection = function( xBag, yBag, xTestBag, 
                             classify, 
                             maxInteractionOrder,
                             nCandidateCovariates,  
                             corFncForCandidateCovariates, corOptionsForCandidateCovariates, 
                             NmandatoryCovariates,
                             interactionsMandatory,
                             keepModel,
                             interactionSeparatorForCoefNames)
{  
  # remove features with missing value or var=0
  indxMiss = apply(is.na(xBag), 2, sum)>0
  indxVar = apply(xBag, 2, var, na.rm=T)==0
  removeFeatures = indxMiss | indxVar
  xBag = xBag[, !removeFeatures, drop=F]
  xTestBag = xTestBag[, !removeFeatures, drop=F]

  nFeatures = ncol(xBag);

  if (NmandatoryCovariates>0)
  {
    mandatCovars = c(1:NmandatoryCovariates)
    # if removed features include mandatory cov, then NmandatoryCovariates should be decreased
    mandatCovars = mandatCovars[ !removeFeatures[mandatCovars]];
    NmandatoryCovariates = length(mandatCovars);
  } else
    mandatCovars = numeric(0);

  nonMandatCovars = setdiff( c(1:nFeatures), mandatCovars);

  # Add interaction terms:
  # generate the interaction matrix
  interactionMatrix = .interactionMatrix(nFeatures, maxInteractionOrder, originalNames = colnames(xBag),
                                    setColNames = TRUE, featureSeparator = interactionSeparatorForCoefNames);
  # Identify mandatory interactions. If interactions of mandatory covariates are also mandatory, add all
  # interactions where at least one term is mandatory (which will be many, so must be used with
  # caution)
  if (interactionsMandatory)
  {
    mandatoryInteractions = apply( interactionMatrix, 2, function(x) { any(x %in% mandatCovars) } );
  } else 
    mandatoryInteractions = mandatCovars;

  nMandatoryInteractions = length(mandatoryInteractions);
  if (nMandatoryInteractions > nCandidateCovariates)
     stop("Number of mandatory interactions is larger than number of candidate covariates.");

  nInteractions = ncol(interactionMatrix);
  nonMandatInteractions = setdiff( c(1:nInteractions), mandatoryInteractions);

  x.int = .generateInteractions(xBag, maxInteractionOrder, 
                                interactionMatrix = interactionMatrix,
                                setColNames = TRUE);
  xTest.int = .generateInteractions(xTestBag, maxInteractionOrder, 
                                interactionMatrix = interactionMatrix,
                                setColNames = TRUE);


  # calculate feature significance  
  corOptionsForCandidateCovariates$x = x.int
  corOptionsForCandidateCovariates$y = yBag
  absGS = abs(do.call(corFncForCandidateCovariates, corOptionsForCandidateCovariates))
  ## nCandidateCovariates could be smaller than indicated due to missing data in x.
  nCandidateCovariates = min(nCandidateCovariates, ncol(x.int))

  ## get indices of candidate cov
  rank = rank(-absGS[nonMandatInteractions], ties.method="f")
  indx = c(mandatoryInteractions, nonMandatInteractions[rank<=(nCandidateCovariates-nMandatoryInteractions)]);

  x.int = x.int[, indx, drop=F]
  xTest.int = xTest.int[, indx, drop=F]
  absGS = absGS[indx]

  # output candidate cov
  candidateFeatures = interactionMatrix[, indx, drop = FALSE];

  # index of most significant feature, used in initial model.
  featureMax = which.max(absGS)

  ## define initial model and full model for binary and continuous outcome. Mandatory covariates must show
  # up in final model, so put them in initial model.

  # Need to remove large variables from the current environment because a copy of this environment is kept
  # in the formula and models below. 
  modelData = data.frame(x.int);

  rm(xBag, xTestBag, x.int, corOptionsForCandidateCovariates, interactionMatrix);

  initialFormula = paste("yBag~", 
               paste(colnames(modelData)[ if (nMandatoryInteractions>0) mandatoryInteractions else featureMax ], 
                     collapse = " + "));
  if (classify) {
    lmInit = glm(formula(initialFormula), data = modelData, family=binomial(link='logit'));
    lmUpper = glm(yBag~., data=modelData, family=binomial(link="logit"))
  } else {
    lmInit = lm(formula(initialFormula), data = modelData);
    lmUpper = lm(yBag~., data = modelData)
  }

  ## forward model selection
  model = stepAIC(lmInit, 
		  scope = list(upper = lmUpper), 
		  direction="forward", 
		  trace=FALSE)

  ## output selected feature (by their names) and their coefficients, is there any smarter way to fish out
  # which features are selected into model? 

  sum = summary(model);
  selected = rownames(sum$coef)[-1]
  featuresInForwardRegression = candidateFeatures[, match(selected, colnames(candidateFeatures)), 
                                                    drop = FALSE];

  coefOfForwardRegression = sum$coef[-1,1]
  interceptOfForwardRegression = sum$coef[1,1]

  ## outHat is piHat for binary outcome and yHat for quantitative outcome
  outHat = predict(model, newdata = as.data.frame(xTest.int), type="response")

  out = list(predicted = outHat, 
             candidateFeatures = candidateFeatures, 
             featuresInForwardRegression = featuresInForwardRegression, 
             coefOfForwardRegression = coefOfForwardRegression, 
             interceptOfForwardRegression = interceptOfForwardRegression,
             model = if (keepModel) model else NULL )

  # The environment of this function is kept in the model returned by stepAIC. Thus, delete everything but
  # the output value so the environment doesn't take up too much memory.
  varList = ls(all.names = TRUE);
  rm(list = setdiff( varList, c("out")));
  out;
}


#=====================================================
#
# Handling of multi-threaded calculations
#
#=====================================================

.disableThreads = function()
{
  try( stopCluster(get(".randomGLMparallelCluster", pos = ".GlobalEnv")), silent = TRUE);
}

.enableThreads = function(nThreads, verbose)
{
  .disableThreads()
  if (is.null(nThreads)) nThreads = max(ceiling(detectCores()*3/4), detectCores()-1);
  if (is.na(nThreads)) nThreads = 1;

  if (nThreads < 1) 
  {
    warning("In function randomGLM: 'nThreads' is below 1. Will use serial execution.");
    nThreads = 1;
  }

  if (nThreads > 1)
  {
    if (verbose > 1) .cat.nl("Will use parallel calculation with ", nThreads, " workers.");
    if (.Platform$OS.type=="windows")
    {
      # On Windows: icreate a cluster manually
      # and  export the parent evinronment of this function for randomGLM to work as well, plus
      # the environment of packages gtools and MASS that are needed. 
      cluster = makePSOCKcluster(nThreads, outfile = "");
      assign(".randomGLMparallelCluster", cluster, pos = ".GlobalEnv");
      clusterExport(cluster, varlist = ls(envir = parent.env(environment()), all.names = TRUE), 
                             envir = parent.env(environment()))
      clusterCall(cluster, library, package = "MASS", character.only = TRUE);
      clusterCall(cluster, library, package = "gtools", character.only = TRUE);
      registerDoParallel(cluster);
    } else {
      # On linux, simply register a parallel backend with nThreads workers
      registerDoParallel(nThreads);
    }
  } else {
    if (verbose > 1) .cat.nl("Will use serial calculation with a single worker process.");
    registerDoSEQ();
  }
  nThreads;
}

#=================================================================================================
#
# randomGLM
#
#=================================================================================================

################################################
## main user level function

randomGLM = function(
  # Input data
  x, y, xtest = NULL, 

  # Include interactions?
  maxInteractionOrder = 1,

  # Prediction type
  classify = is.factor(y) | length(unique(y)) < 4,

  # Multi-level classification options - only apply to classification with multi-level response
  multiClass.global = TRUE,
  multiClass.pairwise = FALSE,
  multiClass.minObs = 1,
  multiClass.ignoreLevels = NULL,

  # Sampling options
  nBags = 100,
  replace = TRUE,
  sampleWeight=NULL,
  nObsInBag = if (replace) nrow(x) else as.integer(0.632 * nrow(x)),
  nFeaturesInBag = ceiling(ifelse(ncol(x)<=10, ncol(x), 
		ifelse(ncol(x)<=300, (1.0276-0.00276*ncol(x))*ncol(x), ncol(x)/5))),
  minInBagObs = min( max( nrow(x)/2, 5), 2*nrow(x)/3),

  # Individual ensemble member predictor options
  nCandidateCovariates=50,
  corFncForCandidateCovariates= cor,
  corOptionsForCandidateCovariates = list(method = "pearson", use="p"),
  mandatoryCovariates = NULL,
  interactionsMandatory = FALSE,
  keepModels = is.null(xtest),

  # Miscellaneous options
  thresholdClassProb = 0.5,
  interactionSeparatorForCoefNames = ".times.",
  randomSeed = 12345,
  nThreads = NULL,
  verbose =0 )
{

  # save original data 
  ySaved = y;
  xSaved = x;
  
  # if y is binary, extract y levels
  if (classify)
  {
    originalYLevels = sort(unique(y));

    # If y has more than 2 levels, do classification on binarized variables. 
    if (length(originalYLevels)>2) 
    {
      if (is.na(multiClass.minObs)) multiClass.minObs = 0;
      if (multiClass.minObs < 1) 
      {
         .cat.nl("Warning: invalid input of 'multiClass.nimObs' changed to 1.");
         multiClass.minObs = 1;
      }
      if (length(originalYLevels) > length(y)/multiClass.minObs | length(originalYLevels) == length(y))
      {
        stop("The response 'y' has too many levels for classification.\n", 
             "   Perhaps you should set 'classify = FALSE'?");
      } else {
        .cat.nl("randomGLM: transforming multi-level response to a series of binary variables.");
      }
   
      yBin = .binarizeCategoricalVar(as.character(y),
                                  minCount = multiClass.minObs, 
                                  val1 = 0, val2 = 1, nameSep = ".vs.", namePrefix = "",
                                  ignore = multiClass.ignoreLevels, 
                                  includePairwise = multiClass.pairwise,
                                  includeLevelVsAll = multiClass.global, 
                                  levelOrder = NULL);
      nY = ncol(yBin);
      yBinNames = colnames(yBin);
      yBinLevels = attr(yBin, "includedLevels");
      
      # Apply randomGLM recursively to each column of yBin.

      out = list(binaryPredictors = list());
      for (iy in 1:nY)
      {
        if (verbose > 0)
          .cat.nl("..Working on binary variable ", yBinNames[iy], " (", iy, " of ", nY, ")");
        out$binaryPredictors[[iy]] = randomGLM(x = x, y = yBin[, iy], xtest = xtest,
                              maxInteractionOrder = maxInteractionOrder,
                              classify = classify,
                              nBags = nBags,
                              replace = replace,
                              sampleWeight = sampleWeight,
                              nObsInBag = nObsInBag,
                              nFeaturesInBag = nFeaturesInBag,
                              nCandidateCovariates = nCandidateCovariates,
                              corFncForCandidateCovariates = corFncForCandidateCovariates,
                              corOptionsForCandidateCovariates = corOptionsForCandidateCovariates,
                              mandatoryCovariates = mandatoryCovariates,
                              interactionsMandatory = interactionsMandatory,
                              keepModels = keepModels,
                              thresholdClassProb = thresholdClassProb,
                              interactionSeparatorForCoefNames = interactionSeparatorForCoefNames,
                              randomSeed = randomSeed,
                              nThreads = nThreads,
                              verbose = verbose - 1);
      }

      names(out$binaryPredictors) = yBinNames;

      out$predictedOOB = as.matrix(sapply(out$binaryPredictors, getElement, "predictedOOB"));
      colnames(out$predictedOOB) = .spaste("PredictionFor.", yBinNames);

      
      out$predictedOOB.response = do.call(cbind, 
                         lapply(out$binaryPredictors, getElement, "predictedOOB.response"));
      responseNames = .spaste(rep(yBinNames, rep(2, nY)), ".ProbabilityOfClass.", as.vector(yBinLevels));
      colnames(out$predictedOOB.response) = responseNames;

      rownames(out$predictedOOB.response) = rownames(out$predictedOOB) = rownames(x);
      if (!is.null(xtest))
      {
        out$predictedTest = sapply(out$binaryPredictors, getElement, "predictedTest");
        colnames(out$predictedTest) = yBinNames;

        out$predictedTest.response = do.call(cbind, 
                         lapply(out$binaryPredictors, getElement, "predictedTest.response"));
        colnames(out$predictedOOB.response) = responseNames;
        rownames(out$predictedTest) = rownames(out$predictedTest.response) = rownames(xtest);
      }

      out$levelMatrix = yBinLevels;
      out$thresholdClassProb = thresholdClassProb;

      class(out) = c("randomGLM", class(out));
     
      return(out);
    }

    y = as.numeric(as.factor(y))-1;
    # The next 3 lines are not needed since the results are always c(0, 1), 0, and 1.
    numYLevels = sort(unique(y));
    minY = min(y, na.rm = TRUE);
    maxY = max(y, na.rm = TRUE);
  } else{
    if (!is.numeric(y)) stop("For quantitative prediction, the response must be numeric.")
  }

  x = as.matrix(x);

  featureNames.original = colnames(x);
  nVars = ncol(x);
  if (is.null(featureNames.original)) 
  {
     featureNames.original = featureNames = .spaste("F", .prependZeros(c(1:nVars)));
     colnames(x) = featureNames.original;
     namesChanged = FALSE;
  } else {
     featureNames = make.names(featureNames.original, unique = TRUE);
     if (isTRUE(all.equal(featureNames, featureNames.original)))
     {
       namesChanged = FALSE;
     } else {
       namesChanged = TRUE;
       nameTranslationTable = data.frame(Column = c(1:nVars), OriginalName = featureNames.original,
                                    CoefficientName = featureNames);
       colnames(x) = featureNames;
     }
  }

  # PL: Why this line? Because this drops all col- and row-names?
  # x = matrix(as.numeric(x), nrow(x), ncol(x))

  if (length(y)!=nrow(x)) {stop("x and y must have the same number of observations.")}

  nSamples = length(y);

  if (nSamples < 8) 
  {
    .cat.nl("*****************************************************************\n",
            "* Warning in randomGLM: there are 7 or fewer observations.\n",
            "*   This may be too few to perform meaningful model selection\n", 
            "*   on in-bag (i.e., even fewer) samples.\n",
            "*   Model selection algorithm will likely output additional warnings.\n",
            "*   The resulting predictor should be used with caution.\n",
            "*****************************************************************");
  }

  nonMandatCovars = setdiff( c(1:nVars), mandatoryCovariates);
  nMandatoryCovariates = length(mandatoryCovariates);


  # nFeaturesInBag shouldn't be greater than total number of features.
  if (nVars<nFeaturesInBag)
  {
    nFeaturesInBag = ncol(x)
    .cat.nl("Warning in randomGLM: nFeaturesInBag is larger than number of features (ncol(x)).\n",
            "   Will use nFeaturesInBag equal to the number of features."); 
  }

  # nCandidateCovariates shouldn't be greater than nFeaturesInBag
  if (nCandidateCovariates>nFeaturesInBag)
  {
     .cat.nl("Warning in randomGLM: nCandidateCovariates is larger than nFeaturesInBag.\n",
             "  Will use nCandidateCovariates=nFeaturesInBag");
     nCandidateCovariates = nFeaturesInBag;
  }

  mandatoryCovarsGiven = !is.null(mandatoryCovariates)
  if (mandatoryCovarsGiven  & nCandidateCovariates < length(mandatoryCovariates))
  {
    stop("Error: number of mandatoryCovariates >= nCandidateCovariates")
  }

  if (thresholdClassProb<0 | thresholdClassProb>1)
    stop("Error: thresholdClassProb takes values between 0  and 1.")

  ## rename features 
  # colnames(x) = paste("feature", 1:nVars, sep="")

  ## scale x
  xSD = apply(x, 2, sd, na.rm=TRUE)
  # xMean = colMeans(x, na.rm=TRUE)
  # x = scale(x)
  validFeatures =xSD>0
  x[, !validFeatures] = 0


  doTest = !is.null(xtest);
  if (doTest)
  { 
    xtestSaved = xtest;

    xtest = as.matrix(xtest);
    if (ncol(x)!=ncol(xtest))
      stop("Number of learning and testing predictors (columns of x, xtest) must equal.");

    if (!is.null(colnames(xtest)))
    {
      if (!isTRUE(all.equal(colnames(xtest), featureNames.original)))
        stop("Column names of 'x' and 'xtest' disagree.");
    } 
    colnames(xtest) = colnames(x)
    nTestSamples = nrow(xtest);
    xtest[, !validFeatures] = 0
   # matrix for test set predicted values across bags 
    predictedTestMat =  matrix(NA, nTestSamples, nBags);
    
  } 

  # set seed
  if (!is.null(randomSeed))
  {
    if (exists(".Random.seed"))
    {
      saved.seed = .Random.seed;
      seedSaved = TRUE;
    } else
      seedSaved = FALSE;
    set.seed(randomSeed);
  }

  # matrix for predicted values in the training set 
  #predictedMat = matrix(NA, nSamples, nBags);
  # matrix for other outputs
  #featuresInForwardRegression = candidateFeatures = coefOfForwardRegression = list();
  #interceptOfForwardRegression = rep(NA, nBags)
  #bagObsIndx = matrix(NA, nBags, nObsInBag)
  #models = list();

  on.exit(.disableThreads())
  nThreads = .enableThreads(nThreads, verbose);

  combinePredictors = function(...)
  {
    preds = list(...);
    out = list();

    out$predictedMat = sapply(preds, getElement, "predicted");
    if (doTest) 
      out$predictedTestMat = sapply(preds, getElement, "predictedTest");

    out$candidateFeatures = lapply(preds, getElement, "candidateFeatures");
    out$featuresInForwardRegression = lapply(preds, getElement, "featuresInForwardRegression");
    out$coefOfForwardRegression = lapply(preds, getElement, "coefOfForwardRegression");
    out$interceptOfForwardRegression = sapply(preds, getElement, "interceptOfForwardRegression");
    out$models = lapply(preds, getElement, "model");
    out;
  }


  # Prepare out of bag samples and in-bag features. This obviates the problems associated with splitting
  # random number generation across worker processes.
  bagFeatures = matrix(NA, nFeaturesInBag, nBags);
  bagObsIndx = matrix(NA, nObsInBag, nBags)
  for (bag in 1:nBags)
  {
    yBagVar = 0;
    while (yBagVar==0)
    {
      # sample indices for each bag
      bagSamples = sample(nSamples, nObsInBag, replace = replace, prob=sampleWeight);
      yBag = y[bagSamples];
      yBagVar = var(yBag, na.rm=TRUE)
      # If there are no out-of-bag samples, force re-sampling as well
      # If the number of in-bag samples is less minInBagObs, re-sample again
      nUniqueInBag = length(unique(bagSamples));
      if (nUniqueInBag==nSamples | nUniqueInBag < minInBagObs) yBagVar = 0
    }
    bagObsIndx[, bag] = bagSamples;
    bagFeatures[, bag] = c(mandatoryCovariates,
                  sample((1:nVars)[nonMandatCovars], nFeaturesInBag - nMandatoryCovariates))
  }

  # For convenience: the iteration over each bag is put in a macro-like function
  singleBagIteration = function(bag, verbose)
  {
    if (verbose>0) {.cat.nl("..bag ", bag)}
    #mem.last = NA;
    #if (verbose > 5) {m = gc()[2,2]; .cat.nl("  Step 1: ", m, ", diff: ", m - mem.last); mem.last = m}
    out = list();
    bagSamples = bagObsIndx[, bag];
    # keep track of in bag and oob samples
    oob = c(1:nSamples)[-unique(bagSamples)];
    nOOB = length(oob);
    features = bagFeatures[, bag]
    xBag = x[bagSamples, features, drop = FALSE];
    yBag = y[bagSamples];

    if (doTest)
    {
      xTestBag = rbind(x[oob, features, drop = FALSE], xtest[, features, drop = FALSE]);
    } else {
      xTestBag = x[oob, features, drop = FALSE];
    }
    # Here I only need to pass the number of mandatory covariates to function forwardSelection, because
    # they're saved at the beginning of features. 

    pr = .forwardSelection(xBag, yBag, xTestBag, 
                          classify=classify, 
                          maxInteractionOrder = maxInteractionOrder,
                          nCandidateCovariates = nCandidateCovariates,
                          corFncForCandidateCovariates = corFncForCandidateCovariates, 
                          corOptionsForCandidateCovariates = corOptionsForCandidateCovariates, 
                          NmandatoryCovariates = nMandatoryCovariates,
                          interactionsMandatory = interactionsMandatory,
                          keepModel = keepModels,
                          interactionSeparatorForCoefNames = interactionSeparatorForCoefNames);

    #if (verbose > 5) {m = gc()[2,2]; .cat.nl("  Step 3: ", m, ", diff: ", m - mem.last); mem.last = m}
    # get output
    out$predicted = rep(NA, nSamples);
    out$predicted[oob] = pr$predicted[1:nOOB];
    if (doTest) {
      out$predictedTest = pr$predicted[(nOOB+1):(nOOB + nTestSamples)];
    }
    # to extract selected feature indices from feature name "feature?", use substring
    dictionary = rbind(c(0,0), cbind(1:nFeaturesInBag, features));
    out$candidateFeatures = .translate(pr$candidateFeatures, dictionary);
    out$featuresInForwardRegression = .translate(pr$featuresInForwardRegression, dictionary);
    out$coefOfForwardRegression = pr$coefOfForwardRegression;
    out$interceptOfForwardRegression = pr$interceptOfForwardRegression

    if (keepModels) out$model = pr$model;
    # For large problems: this may be necessary to preven exhausting the entire memory of the system.
    rm(xBag, pr, xTestBag, features, bagSamples, oob, yBag); 
    #if (verbose > 5) {m = gc()[2,2]; .cat.nl("  Step 4: ", m, ", diff: ", m - mem.last); mem.last = m}

    # Result value for each iteration.
    out
  }

  # loop over bags. Try two different version of the same code, one for parallel and one for serial
  # execution.
  if (nThreads > 1)
  {
    ensemble = foreach (bag = 1:nBags, .combine = combinePredictors, .multicombine = TRUE, 
                                       .maxcombine = nBags) %dopar%
      singleBagIteration(bag, verbose = verbose)
  } else {
    bagRes = list();
    for (bag in 1:nBags)
    {
      bagRes[[bag]] = singleBagIteration(bag, verbose = verbose);
      #if (verbose > 5) {tmp = gc(); .cat.nl("  In main loop: ", tmp[2,2]);}
      #.cat.nl("Size of information for each bag:");
      #print(object.size(bagRes[[bag]]), units = "auto");
    }
    ensemble = do.call(combinePredictors, bagRes);
  }

  # .disableThreads(); <-- this is now included in on.exit() above.

  featuresInForwardRegression.all = do.call(cbind, ensemble$featuresInForwardRegression);
  timesSelectedByForwardRegression = .countsInInteractionMatrix(featuresInForwardRegression.all, nVars)

  colnames(bagObsIndx) = names(ensemble$candidateFeatures) = 
     names(ensemble$featuresInForwardRegression) = names(ensemble$coefOfForwardRegression) = 
     names(ensemble$interceptOfForwardRegression) = .spaste("Bag",1:nBags)

  ## recover original feature names
  if (!is.null(colnames(xSaved))) {
    colnames(timesSelectedByForwardRegression) = colnames(xSaved)
  }

  # average predictive prob over bags, but if all bags give NA, the average should also give NA.
  predictedOOB.response1 = rowMeans(ensemble$predictedMat, na.rm = TRUE);
  predictedOOB.response1[rowSums(!is.na(ensemble$predictedMat))== 0] = NA;
  
  # recover original sample names
  if (!is.null(rownames(xSaved))) {
    names(predictedOOB.response1) = rownames(xSaved)
  }

  # prepare basic output
  out = list(predictedOOB.response = predictedOOB.response1,
             predictedOOB = predictedOOB.response1,
             candidateFeatures= ensemble$candidateFeatures,
             featuresInForwardRegression = ensemble$featuresInForwardRegression,
             coefOfForwardRegression = ensemble$coefOfForwardRegression,
             interceptOfForwardRegression = ensemble$interceptOfForwardRegression,
             bagObsIndx = bagObsIndx,
             timesSelectedByForwardRegression = timesSelectedByForwardRegression,
             models = if (keepModels) ensemble$models else NULL,
             featureNamesChanged = namesChanged,
             nameTranslationTable = if (namesChanged) nameTranslationTable else NULL,
             classify = classify,
             nFeatures = nVars,
             maxInteractionOrder = maxInteractionOrder,
             yLevels = if (classify) originalYLevels else NULL,
	     x.original = xSaved,
	     y.original = ySaved,
             x = x,
             y = y,
             thresholdClassProb = thresholdClassProb)
  # add test set output
  if (doTest) {
    predictedTest.response1 = rowMeans(ensemble$predictedTestMat, na.rm = TRUE);
    predictedTest.response1[rowSums(!is.na(ensemble$predictedTestMat))== 0] = NA;
    if (!is.null(rownames(xtestSaved))) {
       names(predictedTest.response1) = rownames(xtestSaved)
    }

    out$predictedTest.response = predictedTest.response1
    out$predictedTest = predictedTest.response1
  }

  # add output for binary outcomes
  if (classify)
  {
    predictedOOB = ifelse(predictedOOB.response1>thresholdClassProb, 1, 0)
    predictedOOB = originalYLevels[predictedOOB+1]

    predictedOOB.response = cbind(1-predictedOOB.response1, predictedOOB.response1)
    colnames(predictedOOB.response) = as.character(originalYLevels)

    out$predictedOOB = predictedOOB
    out$predictedOOB.response = predictedOOB.response
	
    if (doTest) {
      predictedTest = ifelse(predictedTest.response1>thresholdClassProb, 1, 0)
      predictedTest = originalYLevels[predictedTest+1]

      predictedTest.response = cbind(1-predictedTest.response1, predictedTest.response1)
      colnames(predictedTest.response) = as.character(originalYLevels)

      out$predictedTest.response = predictedTest.response;
      out$predictedTest = predictedTest;
    }
  } 
  class(out) = c("randomGLM", class(out));
  out
}

#============================================================================
#
# predict.randomGLM
#
#============================================================================

# Internal prediction function. This function will also be called from thin.randomGLM to generate prediction
# from the thinned predictor.

.predict.internal = function(object, newdata, type, thresholdClassProb,
                             returnBothTypes = FALSE)
{
  if (!is.null(newdata))
  {
    nSamples = nrow(newdata)
    newdata.1 = cbind(newdata, rep(1, nSamples))
    colnames(newdata) = make.names(colnames(newdata), unique = TRUE);
  } else {
    nSamples = length(object$y);
    x.1 = cbind(object$x, rep(1, nSamples));
  }

  nBags = length(object$models)

  predictedMat = matrix(NA, nSamples, nBags)

  for (b in 1:nBags) if (inherits(object$models[[b]], "lm"))
  {
    bagIM = object$featuresInForwardRegression[[b]];
    if (!is.null(newdata))
    {
      bagNewData = .generateInteractions(x = newdata, x1 = newdata.1, interactionMatrix = bagIM,
                                         maxOrder = object$maxInteractionOrder,
                                         setColNames = TRUE)
      predictedMat[, b] = predict(object$models[[b]], newdata = as.data.frame(bagNewData), 
                                  type = "response");
    } else {
      oob = c(1:nSamples)[-unique(object$bagObsIndx[, b])];
      bagNewData = .generateInteractions(x = object$x[oob, ], x1 = x.1[oob, ], interactionMatrix = bagIM,
                                         maxOrder = object$maxInteractionOrder,
                                         setColNames = TRUE)
      predictedMat[oob, b] = predict(object$models[[b]], newdata = as.data.frame(bagNewData), 
                                  type = "response");
    }
  }

  predicted.response = rowMeans(predictedMat, na.rm = TRUE)
  predicted.response[rowSums(!is.na(predictedMat))== 0] = NA

  names(predicted.response) = if (is.null(newdata)) {
              if (is.null(names(object$y.original))) rownames(object$x.original) else names(object$y.original) 
                                 } else rownames(newdata)

  if (type=="response" | returnBothTypes)
  {
    if (object$classify)
    {
      out.response = cbind(1-predicted.response, predicted.response)
      colnames(out.response) = as.character(object$yLevels)
    }else {
      out.response = predicted.response
    }
  }

  if (type=="class" | returnBothTypes)
  {
    # Note: type == "class" only makes sense for classification. 
    # For continuous prediction, put the continuous prediction here as well.
    if (object$classify)
    {
      predicted.round = ifelse(predicted.response>thresholdClassProb, 1, 0)
      out.class = object$yLevels[predicted.round+1]
    } else
      out.class = predicted.response;
  }

  if (returnBothTypes)
  {
    return(list(response = out.response, class = out.class));
  } else if (type=="class")
    return(out.class);

  out.response;
}
  
#===============================================================================================
#
# user-level predict() function
#
#===============================================================================================
  

predict.randomGLM = function(object, newdata, type=c("response", "class"), 
                             thresholdClassProb = object$thresholdClassProb, ...)
{
  type = match.arg(type)

  if (!is.null(object$binaryPredictors))
  {
    predictions = do.call(cbind, lapply(object$binaryPredictors, predict.randomGLM, 
                            newdata = newdata, type = type, thresholdClassProb = thresholdClassProb, ...))

    if (type=="response")
    {
       colnames(predictions) = colnames(object$predictedOOB.response);
    } else 
       colnames(predictions) = colnames(object$predictedOOB);
    return(predictions);
  }
    
  if (is.null(object$models))
    stop("The 'object' object must contain the undelying models for prediction.\n",
         "   Please re-run the randomGLM function with argument 'keepModels = TRUE' and try again.");

  # If new data is not given, return already calculated prediction. We would need the out-of-bag data for a
  # re-prediction and we don't have them.

  if (missing(newdata))
  {
    stop("valid 'newdata' must be given.")
  }
    
  if (ncol(newdata)!=object$nFeatures)
    stop("Number of columns in 'newdata' differs from the number of features\n",
         "     in the original training data.");


  if (type=="class" & !object$classify)
    stop("type='class' is only valid in classification.")

  if (thresholdClassProb<0 | thresholdClassProb>1)
    stop("Error: thresholdClassProb takes values between 0  and 1.")

  # Call the internal prediction function and return its result.
  .predict.internal(object = object, newdata = newdata, type = type, 
                    thresholdClassProb = thresholdClassProb, returnBothTypes = FALSE);
}

#============================================================================
#
# thin.randomGLM
#
#============================================================================

thinRandomGLM = function(rGLM, threshold)
{

  # Check if rGLM corresponds to multi-level repsonse.
  if (!is.null(rGLM$binaryPredictors))
  {
    out = rGLM;
    out$binaryPredictors = lapply(rGLM$binaryPredictors, thinRandomGLM, 
                               threshold = threshold);

    # Create the main prediction
    yBinNames = names(rGLM$binaryPredictors);
    yBinLevels = rGLM$levelMatrix
    nY = length(rGLM$binaryPredictors);

    names(out$binaryPredictors) = yBinNames;

    out$predictedOOB = as.matrix(sapply(out$binaryPredictors, getElement, "predictedOOB"));
    colnames(out$predictedOOB) = .spaste("PredictionFor.", yBinNames);

    out$predictedOOB.response = do.call(cbind,
                       lapply(out$binaryPredictors, getElement, "predictedOOB.response"));
    responseNames = .spaste(rep(yBinNames, rep(2, nY)), ".ProbabilityOfClass.", as.vector(yBinLevels));
    colnames(out$predictedOOB.response) = responseNames;

    rownames(out$predictedOOB.response) = rownames(out$predictedOOB) = rownames(rGLM$x.original);

    return(out);
  }

  nBags = length(rGLM$featuresInForwardRegression)
  if (nBags != length(rGLM$models))
  if (threshold<0)
    stop("'threshold' must be positive.")

  x = rGLM$x
  y = rGLM$y

  times = rGLM$timesSelectedByForwardRegression[1, ];
  if (threshold >= max(times))
    stop("Specified threshold removes all features in predictor.")

  ## so far, threshold only applies to no interaction level.
  keepF = which(times > threshold)
  screenF = function(input, keepF) { all(is.element(input, keepF)) }

  keepS = rGLM$bagObsIndx
  
  models = featuresInForwardRegression = coefOfForwardRegression = list();
  interceptOfForwardRegression = rep(NA, nBags);
  for (b in 1:nBags)
  {
    yBag = y[ keepS[,b]]
    xBag = x[ keepS[,b],]

    bagIM = rGLM$featuresInForwardRegression[[b]]
    keepIM = apply(bagIM, 2, screenF, keepF)
    if (sum(keepIM)==0)
    {
      featuresInForwardRegression[[b]] = models[[b]] = coefOfForwardRegression[[b]] = NA;
    } else {
      bagIM = bagIM[, keepIM, drop=F]
      featuresInForwardRegression[[b]] = bagIM

      xBag.int = .generateInteractions(x = xBag, interactionMatrix = bagIM, 
	               maxOrder = rGLM$maxInteractionOrder, setColNames = TRUE)

      if (rGLM$classify)
      {
        models[[b]] = glm(yBag~., data=as.data.frame(xBag.int), family=binomial(link="logit"))
      }else{
        models[[b]] = lm(yBag~., data=as.data.frame(xBag.int))
      }
      sum = summary(models[[b]]);
      coefOfForwardRegression[[b]] = sum$coef[-1,1]
      interceptOfForwardRegression[b] = sum$coef[1,1]
    }
  }
  
  featuresInForwardRegression.all = do.call(cbind, featuresInForwardRegression)
  timesSelectedByForwardRegression = .countsInInteractionMatrix(featuresInForwardRegression.all, 
                                                                rGLM$nFeatures)
  if (!is.null(colnames(x))) {
    colnames(timesSelectedByForwardRegression) = colnames(x)
  }

  names(models) = names(featuresInForwardRegression) = names(coefOfForwardRegression) = .spaste("Bag",1:nBags)
  # Copy most of the information from the original rGLM since that information is left intact.
  out = rGLM;
  out$models = models;
  out$featuresInForwardRegression = featuresInForwardRegression;
  out$timesSelectedByForwardRegression = timesSelectedByForwardRegression;
  out$coefOfForwardRegression = coefOfForwardRegression;
  out$interceptOfForwardRegression = interceptOfForwardRegression;

  # Get predictions from the new models, both response and class.

  prediction = .predict.internal(out, newdata = NULL, type = "response",
                                 thresholdClassProb = rGLM$thresholdClassProb,
                                 returnBothTypes = TRUE);


  # Change the elements that need to be changed.
  
  out$predictedOOB.response = prediction$response;
  out$predictedOOB = prediction$class;

  class(out) = c("randomGLM", class(out))
  out
}


