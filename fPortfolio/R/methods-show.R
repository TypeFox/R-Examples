

################################################################################
# FUNCTION:                     DESCRIPTION:
#  show.fPORTFOLIO               S4 Print method for 'fPPORTFOLIO' objects
#  show.fPFOLIODATA              S4 Print method for 'fPFOLIODATA' objects
#  show.fPFOLIOSPEC              S4 Print method for 'fPFOLIOSPEC' objects
#  show.fPFOLIOCON               S4 Print method for 'fPFOLIOCON' objects
################################################################################


setMethod("show", "fPORTFOLIO",
    function(object)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi
    
    # Description:
    #   S4 Print Method for an object of class "fPORTFOLIO"

    # Arguments:
    #   object - an object of class "fPORTFOLIO"

    # FUNCTION:

    # Determine Length Out:
    nFrontierPoints <- NROW(matrix(getWeights(object@portfolio), 
        ncol = getNAssets(object)))
    length.out <- getRmetricsOptions("length.print") # from Rmetrics Options
    index <-
        if (length.out) {
            unique(trunc(seq.int(from = 1, to = nFrontierPoints,
                length.out = length.out)))
        } else {
            seq.int(from = 1, to = NROW(nFrontierPoints))
        }
        
    # Print Title:
    cat("\nTitle:\n ")
        # cat(getType(object),  getTitle(object),     "\n")
        cat(getType(object),  object@title,     "\n")
        cat(" Estimator:        ", getEstimator(object), "\n")
        cat(" Solver:           ", getSolver(object),    "\n")
        cat(" Optimize:         ", getOptimize(object),  "\n")
        cat(" Constraints:      ", getConstraintsTypes(object), "\n")
    if (object@spec@ampl$ampl) {
        cat(" AMPL Project:     ", object@spec@ampl$project, "\n")
        cat(" AMPL Solver:      ", object@spec@ampl$solver, "\n")
    }
    if (!identical(index, 1))
        cat(" Portfolio Points: ", length(index), "of", nFrontierPoints, "\n")
    if (getType(object) == "CVaR")
        cat(" VaR Alpha:        ", getAlpha(object),     "\n")
        #at(" Objective:        ", getObjective(object), "\n")

    # Assets:
    nAssets <- getNAssets(object)
    Names <- names(object@data@statistics$mu)

    # Print Target Weights:
    cat("\nPortfolio Weights:\n")
    table <- matrix(round(getWeights(object@portfolio), digits = 4),
        ncol = nAssets)
    colnames(table) = Names
    rownames(table) = 1:NROW(table)
    print.table(table[index, ])

    # Print Covariance Risk Budgets:
    cat("\nCovariance Risk Budgets:\n")
    table = matrix(round(getCovRiskBudgets(object@portfolio), digits = 4),
        ncol = nAssets)
    colnames(table) = Names
    rownames(table) = 1:NROW(table)
    print.table(table[index, ])

    # PrintCVaR Risk Budgets:
    # to do ...
    
    # Print Tail Risk Budgets:
    # to do ...

    # Print Target Return and Risks:
    # DW: Note object@targetR* is a list do not use getTargetR*()
    targetReturn <- matrix(getTargetReturn(object@portfolio), ncol = 2)
    targetRisk <- matrix(getTargetRisk(object@portfolio), ncol = 4)
    target <- round(cbind(targetReturn, targetRisk), digits = 4)
    if (class(getSeries(object)) == "logical") {
        cat("\nTarget Return and Risk:\n")
        target = target[, c(1, 3), drop = FALSE]
        colnames(target) = c("mean", "Cov")
    } else if( class(getSeries(object)) == "timeSeries") {
        cat("\nTarget Returns and Risks:\n")
        colnames(target) = c("mean", "mu", "Cov", "Sigma", "CVaR", "VaR")
    }
    rownames(target) = 1:NROW(target)
    # DW Only print mu and Sigma for robust estimators!
    if (getEstimator(object) == "covEstimator") {
        print.table(target[index, -c(2,4)])
    } else {
        print.table(target[index, ])
    }
    
    # Print Description:
    cat("\nDescription:\n ")
    # cat(getDescription(object), "\n")
    cat(object@description, "\n")
    
    # Return Value:
    invisible(NULL)
})


# ------------------------------------------------------------------------------


setMethod("show", "fPFOLIODATA",
    function(object)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi
    
    # Description:
    #   S4 Print Method for an object of class "fPFOLIODATA"

    # Arguments:
    #   object - an object of class "fPFOLIOSPEC"

    # FUNCTION:

    # Series:
    cat("\nHead/Tail Series Data:\n")
    if(is.null(dim(object@data$series))) {   
        cat("  No time series data available.\n")
    } else {
        cat("\n")
        print(head(object@data$series, n = 3))
        print(tail(object@data$series, n = 3))
    }
    
    # Statistics:
    cat("\nStatistics:\n\n")
    if(is.null(dim(object@data$series))) { 
        # Print mean and Cov only ..
        print(object@statistics[1:2])
    } else {
        # Print mean, Cov, estimator, mu and Sigma ...
        print(object@statistics)
    }
    
    # Tailrisk:
    # NYI

    # Return Value:
    invisible(NULL)
})


# ------------------------------------------------------------------------------


setMethod("show", "fPFOLIOSPEC",
    function(object)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi
    
    # Description:
    #   S4 Print Method for an object of class "fPFOLIOSPEC"

    # Arguments:
    #   object - an object of class "fPFOLIOSPEC"

    # FUNCTION:

    # Model:

    cat("\nModel List:\t")

    cat("\n Type:                     ",
        object@model$type)

    cat("\n Optimize:                 ",
        object@model$optimize)

    cat("\n Estimator:                ",
        object@model$estimator)

    if (length(object@model$tailRisk) > 0) {
        cat("\n Tail Risk:                ",
            object@model$tailRisk)
    }
    # DW:
    # } else {
    #     cat("\n Tail Risk:                ", "list()")
    # }

    cat("\n Params:                   ",
        paste(names(unlist(object@model$params)), "=",
            unlist(object@model$params)))

    # Portfolio:

    cat("\n\nPortfolio List:\t")

    if (!is.null(object@portfolio$weights)) {
        cat("\n Portfolio Weights:        ",
            object@portfolio$weights)
    } else {
        cat("\n Target Weights:           ", "NULL")
    }

    if (!is.null(object@portfolio$targetReturn)) {
        cat("\n Target Return:            ",
            object@portfolio$targetReturn)
    } else {
        cat("\n Target Return:            ", "NULL")
    }

    if (!is.null(object@portfolio$targetRisk)) {
        cat("\n Target Risk:              ",
            object@portfolio$targetRisk)
    } else {
        cat("\n Target Risk:              ", "NULL")
    }

    if (!is.null(object@portfolio$riskFreeRate)) {
        cat("\n Risk-Free Rate:           ",
        as.character(object@portfolio$riskFreeRate))
    }

    if (!is.null(object@portfolio$nFrontierPoints)) {
        cat("\n Number of Frontier Points:",
        as.character(object@portfolio$nFrontierPoints))
    }

    cat("\n Status:                   ",
        as.character(object@portfolio$status))


    # Optimization:

    cat("\n\nOptim List:\t")

    cat("\n Solver:                   ",
        object@optim$solver)

    if (!is.null(object@optim$objective)) {
        cat("\n Objective:                ",
        as.character(object@optim$objective))
    } else {
        cat("\n Objective:                ", "list()" )
    }

    if (substr(object@optim$solver, 1, 14) != "solveRquadprog") object@optim$options$meq <- NULL
    if (length(object@optim$options) > 0) {
        cat("\n Options:                  ",
            paste(names(unlist(object@optim$options)), "=",
                unlist(object@optim$options)))
    }

    if (length(object@optim$control) > 0) {
        cat("\n Control:                  ",
        as.character(object@optim$control))
    }
    # DW
    # } else {
    #     cat("\n Control:                  ", "list()")
    # }


    cat("\n Trace:                    ",
        object@optim$trace)


    # Messages:

    if (object@messages$messages)
    {
        cat("\n\nMessage List:\t")
    
        if (!is.null(object@messages$list)) {
            cat("\n List:                     ",
                object@messages$list)
        } else {
            cat("\n List:                     ", "NULL")
        }
    }

    # AMPL:
    
    if (object@ampl$ampl)
    {
        cat("\n\nAMPL List:\t")
            
        cat("\n Project:                  ",
            object@ampl$project)
    
        cat("\n Solver:                   ",
            object@ampl$solver)
            
        cat("\n Trace:                    ",
            object@ampl$trace)
    }

    cat("\n")


    # Return Value:
    invisible(NULL)
})


# ------------------------------------------------------------------------------


setMethod("show", "fPFOLIOCON",
    function(object)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi
    
    # Description:
    #   S4 Print Method for an object of class "fPFOLIODATA"

    # Arguments:
    #   object - an object of class "fPFOLIOSPEC"

    # FUNCTION:

    # Print Title:
    cat("\nTitle:\n ")
    cat("Portfolio Constraints\n")

    minmaxW = rbind(object@minWConstraints, object@maxWConstraints)
    rownames(minmaxW) = c("Lower", "Upper")
    if (length(minmaxW)) {
        cat("\nLower/Upper Bounds:\n")
        print(minmaxW)
    }

    eqsumW = object@eqsumWConstraints
    if (sum(dim(eqsumW)) > 2) {
        cat("\nEqual Matrix Constraints:\n")
        print(eqsumW)
    }

    minsumW = object@minsumWConstraints
    if (sum(dim(minsumW)) > 2) {
        cat("\nLower Matrix Constraints:\n")
        print(minsumW)
    }

    maxsumW = object@maxsumWConstraints
    if (sum(dim(maxsumW)) > 2) {
        cat("\nUpper Matrix Constraints:\n")
        print(maxsumW)
    }

    minmaxB <- rbind(object@minBConstraints, object@maxBConstraints)
    if (length(minmaxB) > 0 &&
        !(all(object@minBConstraints == -Inf) && all(object@maxBConstraints == 1))) {
        cat("\nLower/Upper Cov Risk Budget Bounds:\n")
        rownames(minmaxB) = c("Lower", "Upper")
        print(minmaxB)
    }

    listF = object@listFConstraints
    minF = object@minFConstraints
    maxF = object@maxFConstraints
    if (length(listF) > 0) {
        cat("\nNon-Linear Function Constraints:\n")
        minmaxF = rbind(minF, maxF)
        colnames(minmaxF) = names(listF)
        rownames(minmaxF) = c("Lower", "Upper")
        print(minmaxF)
    }
    
    nCard = object@nCardConstraints
    if(nCard > 0) {
        minCard = object@minCardConstraints
        maxCard = object@maxCardConstraints
        minmaxCard = rbind(minCard, maxCard)
        #colnames(minmaxCard) = names(listF)
        rownames(minmaxCard) = c("Lower", "Upper")
        cat("\nCardinality Constraints:\n")
        print(minmaxCard)
    }


    # Return Value:
    invisible(NULL)
})


################################################################################

