###############################################################################
#                        Package: mnlogit                                     #
#                                                                             #
# Multinomial logit, maximum likelihood estimation by Newton-Raphson method   #
#                                                                             #
#               Scientific Computing Group, Sentrana Inc.                     #
###############################################################################

###############################################################################
#                 Main user function                                          # 
# Args:  							              #
#   formula     - same format as mlogit. See help(formula)                    #
#   data        - input data (as a data.frame object) in "long" format        #
#   choiceVar   - the data column containing alternative names                #
#   maxiter     - maximum number of Newton-Raphson iterations to run          #
#   ftol        - function tolerance.                                         #
#                 Difference of two consecutive function evaluation           #
#                 Criteria of terminating Newton's Iterative process          #
#   gtol        - gradient norm tolerance.                                    #
#   weights     - an optional vector of positive frequency weights.           #
#   ncores      - number of processors allowed to use                         #
#   na.rm       - if FALSE then stop(), else remove rows of data with NA      #
#   print.level - increase from 0 to progressively print more runing info     #
#   linDepTol   - tolerance with which linear dep among cols is detected      #
#   start       - initial vector of coefficients                              #
#   alt.subset  - subset of alternatives to perform estimation on             #
#   ...         - currently unused                                            #
#                                                                             #
# Output:                                                                     #
#   mnlogit object                                                            #
#                                                                             #
# Note:                                                                       #
#   'ftol', 'gtol' & 'maxiter' specify Newton-Raphson termination criteria    #
###############################################################################
mnlogit <- function(formula, data, choiceVar=NULL, maxiter = 50, ftol = 1e-6,
             gtol = 1e-6, weights = NULL, ncores = 1, na.rm = TRUE, 
             print.level=0, linDepTol = 1e-6, start=NULL, alt.subset=NULL, ...)
{
    startTime <- proc.time()[3]
    initcall <- match.call()    # Store original function call

    # Basic parameter checking
    if (!is.data.frame(data))
       stop("data must be a data.frame in long format or a mlogit.data object")
    if (ncores < 1) {
        ncores <- 1
        warning("Setting ncores equal to: 1")
    }
    if (!is.null(choiceVar) && is.factor(data[[choiceVar]])) {
        warning(paste("Column", choiceVar, "in data will NOT be treated as",
                      "factor, but as character string!"))
    }

    # Get choiceVar if NULL from data (which MUST now be a mlogit.data object)
    if (is.null(choiceVar) && !any(class(data) == "mlogit.data"))
        stop("Arg data MUST be a mlogit.data object when arg choiceVar = NULL")
    if (is.null(choiceVar)) {
        choiceVar <- "_Alt_Indx_"
        data[[choiceVar]] <- attr(data, "index")$alt # attach to 'data'
    }

    # Extract various types of variables from formula
    formula <- parseFormula(formula) 
    response <- attr(formula, "response")     # response variable
    interceptOn <- attr(formula, "Intercept") # check if intercept is in model
    csvChVar <- attr(formula, "csvChCoeff") 
    indspVar <- attr(formula, "indSpVar") 
    csvGenVar <- attr(formula, "csvGenCoeff") 
    covariates <- c(csvChVar, indspVar, csvGenVar)    
    varNames <- attr(formula, "varNames")    

    if (is.null(covariates) && !interceptOn) 
        stop("Error! Predictor variable(s) must be specified")
    if (is.null(response)) 
        stop("Error! Alternative variable must be specified")
    
    # Handle the alt.subset argument if not NULL 
    if (!is.null(alt.subset)) {
        if (sum(unique(data[[choiceVar]]) %in% alt.subset) < 2)
            stop("Error! Atleast 2 alternatives in data must be in alt.subset")
        keepRows <- data[[choiceVar]] %in% alt.subset
        if (sum(keepRows) <= 0)
            stop("Error! No altrnative in 'alt.subset' is in data.")
        data <- data[keepRows, , drop=FALSE]  
    }

    # Determine relevant parameters
    choice.set <- unique(data[[choiceVar]]) # reordered when data is sorted
    K <- length(choice.set) # number of choices
    if (nrow(data) %% K)
        stop("Mismatch between number of rows in data and number of choices.")
    N <- nrow(data)/K       # number of individuals

    # Check if weights is OK 
    if (!is.null(weights) && length(weights) != N)
      stop("Length of 'weights' arg must match number of observations in data.")
    if (!is.null(weights) && !all(weights > 0))
      stop("All entries in 'weights' must be strictly positive.")      
    # Normalize weights
    if (!is.null(weights)) weights <- weights * N / sum(weights)    

    # Work with only the columns appearing in formula
    data <- data[c(varNames, choiceVar)]

    # Handle NA; Find out row numbers with atleast one NA
    na.rows <- c()
    for (col in 1:ncol(data))
        na.rows <- union(na.rows, which(is.na(data[[col]])))
    Ndropped <- 0
    if (length(na.rows) > 0) {
        if (!na.rm)
            stop("NA present in input data.frame with na.rm = FALSE.")
        # Mark rows with NA for deletion
        keepRows <- rep(TRUE, nrow(data))
        keepRows[na.rows] <- FALSE 
        # Starting with 1st, mark rows for deletion in groups of K
        for (i in 1:N) {
            if (!all(keepRows[((i-1)*K + 1):(i*K)]))
                keepRows[((i-1)*K + 1):(i*K)] <- FALSE
        }
        data <- data[keepRows, , drop=FALSE]
        # Drop weights corresponding to dropped rows 
        if (!is.null(weights)) {
            weights <- weights[keepRows[seq(1, N * K, K)]]
        }
        N <- nrow(data)/K
        Ndropped <- (length(keepRows) - sum(keepRows))/K
    }
    if (print.level && Ndropped > 0) 
      cat(paste("Num of dropped observations (due to NA)  =", Ndropped, "\n"))
      
    # Rearrange the input data.frame object
    # Sort according to choices: data for an atlernative should be contiguous
    data <- data[order(data[[choiceVar]]), ]
    choice.set <- unique(data[[choiceVar]])

    # Obtain response vector as a vector of 0,1
    respVec <- data[[attr(formula, "response")]]
    if (is.factor(respVec)) respVec <- droplevels(respVec)
    respVec <- as.numeric(respVec)
    min.respVec <- min(respVec)
    spread <- max(respVec) - min.respVec
    if (spread != 1) {
        stop(paste("Response variable", attr(formula, "response"), 
                   "must be a factor with exactly 2 levels."))
    }
    respVec <- respVec - min.respVec
    freq.choices <- colSums(matrix(respVec, nrow = N, ncol = K))/N
    loFreq <- min(freq.choices)
    loChoice <- choice.set[which(loFreq == freq.choices)]
    names(freq.choices) <- choice.set
    if (loFreq < 1e-7) {
        cat("Frequencies of alternatives in input data:\n")
        print(prop.table(freq.choices), digits = 4)
        stop(paste("Frequency, in response, of choice:", loChoice, "< 1e-7."))
    }

    # Form design matrices 
    formDesignMat <- function(varVec = NULL, includeIntercept = TRUE)
    {
        if (is.null(varVec) && !includeIntercept) return(NULL) 
        fm <- paste(attr(formula, "response"), "~")
        if (!is.null(varVec))
            fm <- paste(fm, paste(varVec, collapse = "+"))
        if (!includeIntercept) fm <- paste(fm, "-1 ")
        else fm <- paste(fm, "+1 ")
        modMat <- model.matrix(as.formula(fm), data)
    } 
    X <- formDesignMat(varVec = attr(formula, "indSpVar"), 
                       includeIntercept = attr(formula, "Intercept"))
    X <- if (!is.null(X)) X[1:N, , drop=FALSE]   # Matrix of ind sp vars
    Y <- formDesignMat(varVec = attr(formula, "csvChCoeff"), 
                       includeIntercept = FALSE)
    Z <- formDesignMat(varVec = attr(formula, "csvGenCoeff"), 
                       includeIntercept = FALSE)

    # Detect bad columns (these are linearly dependent on other columns)
    badColsList <- list("indSpVar"=NULL, "csvChCoeff"=NULL, "csvGenCoeff"=NULL)
    getNullSpaceCols <- function(mat, tol = 1e-7)
    {
        if (is.null(mat)) return(NULL)
        if (ncol(mat)==1) return(NULL)
        qrdecomp <- qr(mat, tol = tol)
        rank <- qrdecomp$rank
        if (rank == ncol(mat)) return(NULL)
        nullSpCols <- qrdecomp$pivot[(rank + 1):ncol(mat)]
        return(nullSpCols)
    }
    badColsList$indSpVar <- getNullSpaceCols(X, tol = linDepTol)    
    for (i in 1:K) {
        init <- (i-1)*N + 1
        fin <- i*N
        badColsList$csvChCoeff <- union(badColsList$csvChCoeff,
            getNullSpaceCols(Y[init:fin, , drop=FALSE], tol = linDepTol))
    }
    badColsList$csvGenCoeff <- getNullSpaceCols(Z, tol = linDepTol)

    # Get names of variables to be dropped from estimation
    badVarsList <- list()
    badVarsList$indSpVar <- colnames(X[, badColsList$indSpVar, drop=FALSE])
    badVarsList$csvChCoeff <- colnames(Y[, badColsList$csvChCoeff, drop=FALSE])
    badVarsList$csvGenCoeff <- colnames(Z[,badColsList$csvGenCoeff,drop=FALSE])
    badCoeffNames <- makeCoeffNames(badVarsList, choice.set)

    # Eliminate linearly dependent columns
    if (!is.null(X))
      X <- X[ , setdiff(1:ncol(X), badColsList$indSpVar), drop=FALSE]
    if (!is.null(Y))
      Y <- Y[ , setdiff(1:ncol(Y), badColsList$csvChCoeff), drop=FALSE]
    if (!is.null(Z))
      Z <- Z[ , setdiff(1:ncol(Z), badColsList$csvGenCoeff), drop=FALSE]
 
    # Get names of variables 
    varNamesList <- list()
    varNamesList$indSpVar <- colnames(X)
    varNamesList$csvChCoeff <- colnames(Y)
    varNamesList$csvGenCoeff <- colnames(Z)
    coeffNames <- makeCoeffNames(varNamesList, choice.set)
    
    # Order the starting coefficients according to 'coeffNames' 
    if (!is.null(start)) start[coeffNames] <- start
    
    # Do the subtraction: Z_ik - Zi0 (for Generic coefficients data)
    ### NOTE: Base choice (with respect to normalization) is fixed here
    ###       Base choice is the FIRST alternative
    baseChoiceName <- choice.set[1]
    if(!is.null(Z)) { 
        for (ch_k in 2:K) {
            Z[((ch_k - 1)*N + 1):(ch_k*N), ] <-
              Z[((ch_k-1)*N+1):(ch_k*N), , drop=FALSE] - Z[1:N, , drop=FALSE]
        }
    }
    # Drop rows for base alternative
    Z <- Z[(N + 1):(K*N), , drop=FALSE]
    respVec <- respVec[(N + 1):(K*N)]
 
    t1 <- proc.time()[3]    # Time at end of pre-processing

    gc()  # Invoke garbage collector at end of pre-processing 
    prep.time <- t1 - startTime
    if (print.level > 1) {
      cat(paste0("Base alternative is: ", baseChoiceName))
      cat(paste0("\nPreprocessing data for estimation took ", 
                  round(prep.time, 3), " sec.\n"))
    } 
    # Solve MLE using Newton-Raphson
    result <- newtonRaphson(respVec, X, Y, Z, K, maxiter, gtol, ftol, ncores,
                  print.level, coeffNames, weights=weights, start=start)
    result$est.stats$prepTimeSecs <- prep.time
    # Post-processing
    colnames(result$hessMat) <- coeffNames 
    rownames(result$hessMat) <- coeffNames
    names(result$grad)   <- coeffNames

    # Reorder coeff so that those for a choice are together 
    od <- reordering(varNamesList, choice.set)
    coeffNames <- makeCoeffNames(varNamesList, choice.set)
    coefficients <- c(result$coeff, if (is.null(badCoeffNames)) NULL
                                    else rep(NA, length(badCoeffNames)))
    names(coefficients) <- c(coeffNames,
        badCoeffNames[reordering(badVarsList, choice.set)])
    reordered_coeff <- c(result$coeff[od], if (is.null(badCoeffNames)) NULL
                                           else rep(NA, length(badCoeffNames)))
    names(reordered_coeff) <- c(coeffNames[od],
        badCoeffNames[reordering(badVarsList, choice.set)])

    # Set colnames for probability & residual matrix
    colnames(result$probability) <- choice.set
    if (maxiter > 0) colnames(result$residual)    <- choice.set

    result$model.size$intercept <- interceptOn
    attributes(formula) <- NULL
   
    # Loglikelihood and AIC  
    logLik <- structure(-result$loglikelihood,
                        df = result$model.size$nparams,    
                        class = "logLik"
                       )
    AIC <- 2*(result$model.size$nparams + result$loglikelihood)

    # 'index' attribute for data
    index <- data.frame(chid = rep(1:result$model.size$N, result$model.size$K),
                        alt = data[[choiceVar]])
    attr(data, "index") <- index
    
    fit <- structure(
             list(
               coefficients  = coefficients,
               logLik        = logLik,
               gradient      = -result$grad,
               hessian       = result$hessMat,
               est.stat      = result$est.stats,
               fitted.values = 1 - attr(result$residual, "outcome"),
               probabilities = result$probability,
               residuals     = result$residual,
               df            = result$model.size$nparams, #
               AIC           = AIC,                       #
               choices       = choice.set,                # 
               model.size    = result$model.size,         #
               ordered.coeff = reordered_coeff,           #
               model         = data,                     
               freq          = freq.choices,             
               formula       = Formula(formula(formula)),
               call          = initcall),
             class = "mnlogit"
           )

    if (print.level)
      cat(paste0("\nTotal time spent in mnlogit = ",
                 round(proc.time()[3] - startTime, 3), " seconds.\n"))
    return(fit)
}

# Makes names of model coefficients
# Sets ordering of the Hessian and gradient rows.
# Ensures order is in confirmation with design matrices 
makeCoeffNames <- function (varNames, choices)
{
    if (length(varNames) == 0) return(NULL)
    choices <- as.vector(choices)
    coeffName <- c(outer(varNames$indSpVar, choices[-1], paste, sep=":"),
                   outer(varNames$csvChCoeff, choices, paste, sep=":"),
                   varNames$csvGenCoeff)
}

# Generate a re-ordering of coeff names (group choices together)
reordering <- function(varList, choices)
{
    if (length(varList) == 0) return(NULL)
    K <- length(as.vector(choices))
    p <- length(as.vector(varList$indSpVar))
    f <- length(as.vector(varList$csvChCoeff))
    d <- length(as.vector(varList$csvGenCoeff))
    orig <-  c(if (p > 0) rep(1:p, K-1) else NULL,
               if (f > 0) rep((p+1):(p+f), K) else NULL,
               if (d > 0) (p+f+1):(p+f+d) else NULL)
    order(orig) 
}
