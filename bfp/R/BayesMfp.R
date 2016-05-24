#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[BayesMfp.R] by DSB Mit 30/05/2012 13:32 (CEST)>
##
## Description:
## Main function of the bfp package: Bayesian Inference for a multivariate
## fractional polynomial model via "BayesMfp" function. Plus model formula helper
## functions "bfp" and "uc".
##
## History:
## 22/08/2008   copy from BayesMfp package.
## 23/08/2008   modified according to new methodology:
##              - prior specification is now only the hyperparameter a
## 03/09/2008   - no changes are necessary for new elements in return list
##              - enforce that only the intercept can be a fixed term in the formula
## 04/09/2008   - add SST and yMean attribute to return list
## 05/09/2008   - add xCentered attribute to return list
## 13/11/2008   add rangeVals parameter to bfp, which is useful for following prediction
##              on new covariates data, which would otherwise not be possible if the
##              ranges differ too much from the training data.
## 14/11/2008   save .Random.seed at the beginning,
##              now correct scaling if scale=FALSE for some bfp term
## 29/11/2008   new prior specification option: flat prior on model space or old prior
##              which favours sparse models (more strongly)
## 21/09/2009   Check that only for k > 1 the model sampling approach can be chosen
##              (necessary because of the new SWITCH move type), and that all
##              maximum FP degrees are identical (to be inline with the paper)
## 05/10/2009   remove dead comments
## 03/12/2009   add option to specify the models cache size when doing sampling
## 14/12/2009   increase default value for "nCache" to 1e9 (previously 1e5) for
##              (almost) backward-compatibility - if not more than 1e9 models are visited
##              during the sampling, the results are completely identical to those from
##              the old HypergBayesMfp package.
## 29/10/2010   add additional attribute "linearInclusionProbs", which gives for
##              FP terms the posterior probability of exactly linear inclusion.
## 14/01/2011   - add the new modelPrior option "dependent"
##              - let getNumberPossibleFps be a separate function so that it can
##              be used in the function getLogPrior as well (for the dependent
##              model prior)
#####################################################################################

getNumberPossibleFps <- function (  # computes number of possible univariate fps (including omission)
                                  maxDegree # maximum fp degree
                                  ){
    s <- ifelse (maxDegree <= 3, 8, 5 + maxDegree) # card. of power set
    singleDegreeNumbers <- sapply (0:maxDegree, function (m)
                                       choose (s - 1 + m, m))
    return (sum (singleDegreeNumbers))
}

`BayesMfp` <-
    function (
              formula = formula(data), # model formula
              data = parent.frame(),   # data.frame for model variables
              family = gaussian,       # distribution and link: only gaussian supported at the moment
              priorSpecs =             # prior specifications:
              list (a = 4,             # hyperparameter for hyper-g prior, must be greater than 3
                    modelPrior = "flat"), # prior on model space: "flat",
                                        # "sparse" or "dependent"
              method = c ("ask", "exhaustive", "sampling"), # which method should be used to explore the
                                        # posterior model space? default is to ask after prompting the space cardinality
              subset = NULL,           # optional subset expression
              na.action = na.omit,     # default is to skip rows with missing data
              verbose = TRUE,          # should information on computation progress be given?
              nModels = NULL,          # how many best models should be saved? default: 1% of total
                                       # number of (cached) models. Must not be larger than nCache
                                       # if method == "sampling".
              nCache=1e9L,              # maximum number of best models to be cached at the same
                                        # time during the model sampling, only has effect if method = sampling
              chainlength = 1e5L        # only has effect if method = sampling
              )
{
    ## save call for return object
    call <- match.call()
    method <- match.arg (method)

    ## save random seed
    randomSeed <- 
        if(exists(".Random.seed", .GlobalEnv))
            get(".Random.seed", .GlobalEnv)
        else
            NULL
    
    ## evaluate family
    if (is.character(family))
        family <- get(family, mode = "function")
    if (is.function(family))
        family <- family()
    if (family$family != "gaussian" | family$link != "identity")
        stop (simpleError("currently only the Gaussian distribution with identity link is implemented"))
    if (is.null(family$family)) {
        print(family)
        stop(simpleError("`family' not recognized"))
    }

    ## evaluate call for model frame building
    m <- match.call(expand.dots = FALSE)

    ## select normal parts of the call
    temp <- c("", "formula", "data", "subset", "na.action") # "" is the function name
    m <- m[match(temp, names(m), nomatch = 0)]

    ## sort formula, so that bfp comes before uc
    Terms <- if (missing(data))
        terms(formula)
    else
        terms(formula, data = data)

    ## check if intercept is present
    if (! attr(Terms, "intercept"))
        stop(simpleError("There must be an intercept term in the model formula."))

    ## now sort the formula
    sortedFormula <- paste (deparse (Terms[[2]]),
                            "~ 1 +", 
                            paste (sort (attr (Terms, "term.labels")), collapse = "+")
                            )
    sortedFormula <- as.formula (sortedFormula)

    ## filter special parts in formula: uncertain covariates (uc) and (Bayesian) fractional polynomials (bfp)
    special <- c("uc", "bfp")
    Terms <- if (missing(data))
        terms(sortedFormula, special)
    else
        terms(sortedFormula, special, data = data)

    ucTermInd <- attr (Terms, "specials")$uc  # special indices in original formula (beginning with 1 = response!)
    nUcGroups <- length (ucTermInd)
    bfpTermInd <- attr (Terms, "specials")$bfp
    nFps <- length (bfpTermInd)
    ## check if bfp's are present
    if (! nFps)
        stop (simpleError("No fractional polynomial terms in formula! Aborting."))

    vars <- attr (Terms, "variables")    # language object
    varlist <- eval(vars, envir = data)              # list
    covariates <- paste(as.list (vars)[-c(1,2)]) # vector with covariate entries (no list or response or Intercept)

    bfpInner <- varlist[bfpTermInd]     # saved for later use
    covariates[bfpTermInd - 1] <- unlist(bfpInner) # remove bfp( ) from formula; -1 because of reponse column

    if (nUcGroups){                     # if ucs are present
        ucInner <- unlist(varlist[ucTermInd])
        covariates[ucTermInd - 1] <- ucInner
        ## determine association of terms with uc groups
        ucTermLengths <- sapply (ucInner, function (oneUc)
                                 length (attr (terms (as.formula (paste ("~", oneUc))), "term.labels"))
                                 )
        ucTermLengthsCum <- c(0, cumsum (ucTermLengths - 1)) # how much longer than 1, accumulated
        ucTermList <- lapply (seq (along = ucTermInd), function (i) # list for association uc group and assign index
                              as.integer(
                                         ucTermInd[i] - 1 + # Starting assign index
                                         ucTermLengthsCum[i]+ # add lengths from before
                                         0:(ucTermLengths[i]-1) # range for this uc term
                                         )
                              )
    } else {
        ucInner <- ucTermList <- NULL
    }

    ## build new Formula
    newFormula <-                       # is saved for predict method at the end
        update (sortedFormula,
                paste (".~ 1 +",  
                       paste (covariates, collapse = "+")) # only update RHS
                )
    newTerms <- if (missing(data))
        terms(newFormula)
    else
        terms(newFormula, data = data)

    ## build model frame
    m$formula <- newTerms
    m$scale <- m$family <- m$verbose <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval(m, envir=sys.parent())
    Y <- model.response(m, type = "double")

    ## build design matrix
    X <- model.matrix (newTerms, m)
    Xcentered <- scale(X, center=TRUE, scale=FALSE)

    ## which terms gave rise to which columns?
    ## (0 = intercept, 1 = first term)
    termNumbers <- attr (X, "assign") 

    ## vector of length col (X) giving uc group indices or 0 (no uc)
    ## for associating uc groups with model matrix columns: ucIndices
    ucIndices <- fpMaxs <- integer (length (termNumbers))

    ## list for mapping group -> columns in model matrix: ucColList
    if(nUcGroups){
        for (i in seq (along=ucTermList)){
            ucIndices[termNumbers %in% ucTermList[[i]]] <-  i
        }

        ucColList <- lapply (seq (along=ucTermList), function (ucGroup) which (ucIndices == ucGroup))
    } else {
        ucColList <- NULL
    }


    ## vectors of length col (X) giving maximum fp degrees or 0 (no bfp)
    ## and if scaling is wanted (1) or not (0)
    ## for associating bfps with model matrix columns
    ## In addition, scale Columns or exit if non-positive values occur
    bfpInds <- bfpTermInd - 1 + attr (Terms, "intercept")           # now indexing matrix column
    for (i in seq (along=bfpInner)){
        colInd <- bfpInds[i]
        fpObj <- bfpInner[[i]]

        fpMaxs[colInd] <- attr (fpObj, "max")

        ## get scaling info
        scaleResult <- fpScale (c(attr(fpObj, "rangeVals"), # extra values not in the data
                                  X[,colInd]),              # covariate data
                                scaling = attr(fpObj, "scale")) # scaling wished?
        attr (bfpInner[[i]], "prescalingValues") <- scaleResult
        
        ## do the scaling
        X[,colInd] <- X[,colInd] + scaleResult$shift
        X[,colInd] <- X[,colInd] / scaleResult$scale

        ## check positivity
        if (min(X[,colInd]) <= 0) 
            stop (simpleError(paste("Prescaling necessary for negative values in variable ", fpObj,
                                    "! Aborting.")))
    }

    ## check that all maximum FP degrees are equal, so that the SWITCH will always be possible.
    ## This assumption could potentially be removed later on, or otherwise the max option in bfp
    ## could be removed.
    if (! identical(length(unique(fpMaxs[fpMaxs != 0])),
                    1L))
        stop(simpleError("All maximum FP degrees must be identical."))

    ## check if only the intercept (one column) is a fixed term
    if (sum(! (fpMaxs | ucIndices)) > 1)
        stop (simpleError("Only the intercept can be a fixed term."))

    ## check hyperparameter for hyper-g
    if (priorSpecs$a <= 3){
        stop(simpleError(paste("Hyperparameter for hyper-g prior leads to infinite posterior expected g's:\n",
                               "a must be greater than 3")))
    } else if (priorSpecs$a > 4) {
        warning(simpleWarning(paste("Hyperparameter", priorSpecs$a,
                                    "for hyper-g prior puts more mass",
                                    "on shrinkage values\n",
                                    "near zero than the uniform prior on the shrinkage factor. This may",
                                    "be undesirable.")))
    }

    ## get model prior choice
    priorSpecs$modelPrior <- match.arg(priorSpecs$modelPrior,
                                       choices=c("flat", "sparse", "dependent"))

    ## compute and print cardinality of the model space to guide decision
    fpSetCards <- ifelse (fpMaxs[fpMaxs != 0] <= 3, 8, 5 + fpMaxs[fpMaxs != 0])
    singleNumbers <- sapply (fpMaxs, getNumberPossibleFps)
    totalNumber <- prod (singleNumbers) * 2^(nUcGroups) # maximum number of possible models
        
    ## decide method
    if (identical(method, "ask")){
        cat ("The cardinality of the model space is at most ", totalNumber, ".\n", sep = "")
        decision <- substr (readline(paste("Do you want to do a deterministic search for the best model (y)",
                                           "or sample from the model space (n) or abort (else) ?\n")),
                            1, 1)
    } else {
        decision <- switch (method, exhaustive = "y", sampling = "n")
    }

    ## ensure that we only do model sampling if there is more than 1 FP term in the model
    if(identical(decision, "n") && identical(nFps, 1L))
    {
        warning(simpleWarning(paste("We need to do an exhaustive computation of all models,",
                                    "because there is only 1 FP term in the model!")))
        decision <- "y"
    }

    ## start

    if (identical(decision, "n")){
        ## get chainlength?
        if (identical(method, "ask"))
        {
            chainlength <- as.numeric (readline ("How long do you want the Markov chain to run?\n"))
        }

        ## compute the default number of models to be saved
        if(is.null(nModels))
        {
            nModels <- max(chainlength / 100, 1)
        }
        else
            stopifnot(nModels >= 1)
            
        
        ## check the chosen cache size
        nCache <- as.integer(nCache)
        stopifnot(nCache >= nModels)        

        ## echo progress?
        if (verbose){
            cat("Starting sampler...\n")
            cat ("0%", rep ("_", 100 - 6), "100%\n", sep = "")
        }

        ## then go C++
        Ret <-
            .Call ("samplingGaussian", ## PACKAGE = "bfp",
                   X,            # design matrix
                   Xcentered,    # centered design matrix
                   Y,            # response vector
                   as.integer(fpMaxs[fpMaxs != 0]), # vector of maximum fp degrees
                   as.integer(bfpInds), # vector of fp columns
                   as.integer (fpSetCards), # cardinality of corresponding power sets
                   as.integer(nFps), # number of fp terms
                   unlist (bfpInner), # names of fp terms
                   as.integer(ucIndices), # vector giving uncertainty custer indices (column -> which group)
                   ucColList,   # list for group -> which columns mapping
                   as.integer(nUcGroups), # number of uncertainty groups
                   as.double (priorSpecs$a), # only the hyperparameter a
                   priorSpecs$modelPrior, # model prior?
                   as.integer(nModels),          # number of best models returned
                   verbose,          # should progress been displayed?
                   as.double(chainlength), # how many times should a jump be proposed?
                   as.integer(nCache)      # size of models cache (an STL map)
                   )

        attr (Ret, "chainlength") <- chainlength

    } else if (identical(decision, "y")){

        ## compute the default number of models to be saved
        if(is.null(nModels))
        {
            nModels <- max(totalNumber / 100, 1)
        }
        else
            stopifnot(nModels >= 1)

        ## echo progress?
        if (verbose){
            cat("Starting with computation of every model...\n")
            cat ("0%", rep ("_", 100 - 6), "100%\n", sep = "")
        }

        ## then go C++
        Ret <-
            .Call ("exhaustiveGaussian", ##  PACKAGE = "bfp",
                   X,            # design matrix
                   Xcentered,    # centered design matrix
                   Y,            # response vector
                   as.integer(fpMaxs[fpMaxs != 0]), # vector of maximum fp degrees
                   as.integer(bfpInds), # vector of fp columns
                   as.integer (fpSetCards), # cardinality of corresponding power sets
                   as.integer(nFps), # number of fp terms
                   unlist (bfpInner), # names of fp terms
                   as.integer(ucIndices), # vector giving uncertainty custer indices (column -> which group)
                   ucColList,   # list for group -> which columns mapping
                   as.integer(nUcGroups), # number of uncertainty groups
                   as.double(totalNumber), # cardinality of model space
                   as.double (priorSpecs$a), # only the hyperparameter a
                   priorSpecs$modelPrior, #  model prior?
                   as.integer(nModels),          # number of best models returned
                   verbose          # should progress been displayed?
                   )

    } else {

        cat ("Aborting.\n")
        return ()
    }

    ## C++ attaches the following attributes:

    ## numVisited 
    ## inclusionProbs
    ## linearInclusionProbs (only for FPs)
    ## logNormConst

    names (attr (Ret, "inclusionProbs")) <- c(unlist (bfpInner), ucInner)
    names (attr (Ret, "linearInclusionProbs")) <- c(unlist (bfpInner))

    ## attach additional information:
    names (Ret) <- 1:length (Ret)

    attr (Ret, "call") <- call
    attr (Ret, "formula") <- newFormula

    attr (Ret, "x") <- X
    attr(Ret, "xCentered") <- Xcentered
    attr (Ret, "y") <- Y
    
    attr(Ret, "SST") <- sum((Y - mean(Y))^2)
    attr(Ret, "yMean") <- mean(Y)
    
    fixedInds <- setdiff (1:ncol (X), c (bfpInds, which (ucIndices > 0)))
    attr (Ret, "indices") <- list (uc = ucIndices, ucList = ucColList, bfp = bfpInds, fixed = fixedInds)

    fixedNamesInds <- setdiff(2:length (varlist), unlist (attr (Terms, "specials"))) + 1
    interceptName <- ifelse (attr (Terms, "intercept"), "(Intercept)", NULL)
    attr (Ret, "termNames") <- list (fixed = c(interceptName, sapply (as.list (vars[fixedNamesInds]), deparse)),
                                     bfp = unlist (bfpInner),
                                     uc = ucInner
                                     )

    shiftScaleMaxMat <- matrix (nrow = nFps, ncol = 4)
    shiftScaleMaxMat[, 1:2] <- matrix(unlist(lapply (bfpInner, attr, "prescalingValues")), ncol = 2, byrow = TRUE)
    shiftScaleMaxMat[, 3] <- fpMaxs[fpMaxs != 0]
    shiftScaleMaxMat[, 4] <- fpSetCards

    rownames (shiftScaleMaxMat) <- unlist (bfpInner)
    colnames (shiftScaleMaxMat) <- c ("shift", "scale", "maxDegree", "cardPowerset")

    attr (Ret, "shiftScaleMax") <- shiftScaleMaxMat
    attr (Ret, "priorSpecs") <- priorSpecs
    attr(Ret, "randomSeed") <- randomSeed
    
    ## and return
    class (Ret) <- c("BayesMfp", "list")
    return(Ret)
}

####################################################################################################

`bfp` <-
function (
          x,                    # variable
          max = 2,              # maximum degree for this fp
          scale = TRUE,         # use pre-transformation scaling to avoid numerical problems?
          rangeVals = NULL      # extra numbers if the scaling should consider values in this range
          )
{
    x <- deparse(substitute(x))
    attr(x, "max") <- max
    attr(x, "scale") <- scale
    attr(x, "rangeVals") <- rangeVals
    x
}

####################################################################################################

`uc` <-
function (x)
{
    x <- deparse(substitute(x))
    x
}
