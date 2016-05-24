

################################################################################
# FUNCTION:                DESCRIPTION:
#  portfolioSpec            Returns an object of class fPFOLIOSPEC
#  .checkWeights            Checks and forces tiny weights to zero
#  .checkSpecVsConstraints  Checks if spec and constraints do match
#  .checkTargetReturn       Checks if target Return is defined
################################################################################


portfolioSpec <-
function(
    model = list(
        type = "MV",                   # Alt: "LPM", "CVaR"
        optimize = "minRisk",          # Alt: "maxReturn"
        estimator = "covEstimator",    # Alt: "shrinkEstimator", 
                                       #      "lpmEstimator"
        tailRisk = list(),
        params = list(alpha = 0.05, a = 1)),
    portfolio = list(
        weights = NULL,
        targetReturn = NULL,
        targetRisk = NULL,
        riskFreeRate = 0,
        nFrontierPoints = 50,
        status = NA),
    optim = list(
         solver = "solveRquadprog",     # Alt: "solveRdonlp2" 
                                        #      "solveRsolnp"
                                        #      "solveRglpk", 
                                        #      "solveRsymphony"
                                        #      "solveRsocp" ...
         objective = c(
             "portfolioObjective", 
             "portfolioReturn", 
             "portfolioRisk"),
         options = list(meq = 2),
         control = list(),
         trace = FALSE),
    messages = list(
        messages = FALSE,
        note = ""),
    ampl = list(
        ampl = FALSE,
        project = "ampl",
        solver = "ipopt",
        protocol = FALSE,
        trace = FALSE)
    )
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Specifies a portfolio to be optimized

    # Example:
    #   portfolioSpec(portfolio = list(targetReturn = 1.5))

    # FUNCTION:

    # Compose Checklists:
    # model.type = c("MV", "CVaR")
    # model.estimator.mean = "mean"
    # model.estimator.cov = c("cov", "mcd", "Mcd", "shrink")
    # optim.solver = c("solveRquadprog", "solveRdonlp2", "solveRglpk")
    # optim.trace = FALSE

    # Check Arguments:
    # stopifnot(model$type %in% model.type)
    # stopifnot(model$estimator[1] %in% model.estimator.mean)
    # stopifnot(model$estimator[2] %in% model.estimator.cov)
    # stopifnot(optim$solver %in% optim.solver)

    # Model Slot:
    Model = list(
        type = "MV",
        optimize = "minRisk",
        estimator = "covEstimator",
        tailRisk = NULL,
        params = list())
    model$type = model$type[1]
    Model[(Names <- names(model))] <- model

    # Portfolio Slot:
    Portfolio = list(
        weights = NULL,
        targetReturn = NULL,
        targetRisk = NULL,
        riskFreeRate = 0,
        nFrontierPoints = 50,
        status = 0)
    Portfolio[(Names <- names(portfolio))] <- portfolio

    # Check Portfolio - weights, targetReturn, targetRisk:
    # ... at least two of them must be set to NULL!
    checkPortfolio = 0
    if(!is.null(portfolio$weights)) checkPortfolio = checkPortfolio + 1
    if(!is.null(portfolio$targetReturn)) checkPortfolio = checkPortfolio + 1
    stopifnot(checkPortfolio <= 1)

    # Optim Slot:
    Optim = list(
        solver = "solveRquadprog",
        objective = NULL,
        options = list(meq = 2),
        control = list(),
        trace = FALSE)
    Optim[(Names <- names(optim))] <- optim
    
    # Messages Slot:
    Messages = list(
        list = NULL)
    Messages[(Names <- names(messages))] <- messages

    # Return Value:
    new("fPFOLIOSPEC",
        model = Model,
        portfolio = Portfolio,
        optim = Optim,
        messages = messages,
        ampl = ampl)
}


# ------------------------------------------------------------------------------


.checkSpec <-
function(spec)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Checks for specification conflicts
    
    # FUNCTION:
    
    if (getSolver(spec) == "solveRglpk" && getType(spec) == "MV") {
        # Error Message:
        cat("\nExecution stopped:")   
        cat("\n  Specification conflict for portfolio solver and type.")
        cat("\nSpec Information:")
        cat("\n  Solver=", getSolver(spec), ",", " type = ", getType(spec), 
            ".", sep = "")
        cat("\n")
        stop(call. = FALSE, show.error.messages = "\n  returned from Rmetrics")
     }
     
     if (getSolver(spec) == "solveRsymphony" && getType(spec) == "MV") {
        # Error Message:
        cat("\nExecution stopped:")   
        cat("\n  Specification conflict for portfolio solver and type.")
        cat("\nSpec Information:")
        cat("\n  Solver=", getSolver(spec), ",", "type = ", getType(spec), 
            ".", sep = "")
        cat("\n")
        stop(call. = FALSE, show.error.messages = "\n  returned from Rmetrics")
     }
     
     # Return Value:
     "ok"
}


# ------------------------------------------------------------------------------


.checkWeights <-
function(weights, eps = sqrt(.Machine$double.eps))
{    
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Sets tiny weights to zero
    
    # Arguments:
    #   weights - a numeric vector of portfolio weights
    #   eps - a numeric value, lower bounds of weigths
    
    # FUNCTOION:
    
    # Check:
    for(i in 1:length(weights)) {
        if(abs(weights[i]) < eps) weights[i] = 0
    }
    
    # Return Value:
    weights
}


# ------------------------------------------------------------------------------


.checkSpecVsConstraints <-
function(spec, constraints)
{    
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Check if spec and constraints do match
    
    # Arguments:
    #   spec - portfolio specification as fPFOLIOSPEC object
    #   constraints - as charvec or as fPFOLIOSPEC object
    
    # FUNCTOION:
    
    # Check:
    if(class(constraints) == "fPFOLIOCON")
        constraints = constraints@stringConstraints
    if(any(constraints == "Short")) {
        stopifnot(getSolver(spec) == "solveRshortExact")
    }
    
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


.checkTargetReturn <-
function(spec)
{    
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Check if target Return is defined
    
    # Arguments:
    #   spec - specification object
    
    # FUNCTOION:
    
    # Check:
    targetReturn = getTargetReturn(spec)
    if(is.null(targetReturn))
        stop("The target return is not available")
    
    # Return Value:
    invisible(targetReturn)
}


################################################################################

