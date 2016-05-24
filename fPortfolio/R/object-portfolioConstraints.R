

################################################################################
# FUNCTION:               DESCRIPTION:
#  portfolioConstraints    Returns an object of class fPFOLIOCON
# FUNCTION:               DESCRIPTION:
#  minWConstraints         Returns vector with min box constraints
#  maxWConstraints         Returns vector with max box constraints
#  eqsumWConstraints       Returns list with group equal vec/matrix constraints
#  minsumWConstraints      Returns list with group min vec/matrix constraints
#  maxsumWConstraints      Returns list with group max vec/matrix constraints
#  minBConstraints         Returns vector with min cov risk budget constraints
#  maxBConstraints         Returns vector with max cov risk budget constraints
#  minFConstraints         Returns vector with min nonlin functions constraints
#  maxFConstraints         Returns vector with max nonlin functions constraints
#  minBuyinConstraints     Returns lower bound of buyin constraints
#  maxBuyinConstraints     Returns upper bound of buyin constraints
#  nCardConstraints        Returns number of Cardinalities
#  minCardConstraints      Returns lower bound of Cardinalities
#  maxCardConstraints      Returns upper bound of Cardinalities
################################################################################


portfolioConstraints <-
function(data, spec=portfolioSpec(), constraints="LongOnly", ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns an object of class fPFOLIOCON

    # Arguments:
    #   data - a timeSeries or a fPFOLIODATA object
    #   spec - a fPFOLIOSPEC object
    #   constraints - a constraints string
    #       validStrings = c(
    #           "LongOnly", "Short",      # LongOnly and Short Notification
    #           "minW", "maxW",           # Box Constraints
    #           "minsumW", "maxsumW",     # left/right Sided Group Constraints
    #           "minB", "maxB",           # Covariance Risk Budgets
    #           "listF", "minF", "maxF",  # Nonlinear Functions Constraints
    #     NEW:  "minBuyin"," maxBuyin",
    #     NEW:  "nCard", "minCard", "maxCard") 
    
    # Details:
    #   This function takes the constraints strings and converts them to
    #   constraints vectors and matrices of the following form:
    #   1. boxConstraints           W_min <= W <= W_max
    #   2. groupEqConstraints       A_eq W = c_eq
    #   3. groupMatConstraints      a_vec <= A_mat W <= b_vec
    #   4. riskBudgetConstraints    a <= RiskBudget <= b
    #      cardinalityConstraints   eps*z <= W <- delta*z,  z[0,1],  Sum(z)=K
    #   These values are returned as list in four slots.

    # Example:
    #   data = .lppData; spec=.mvSpec
    #   portfolioConstraints(data, spec, "LongOnly")
    #   constraints=c("minW[1:3]=0.1", "maxW[4:6]=0.9", "minsumW[c(2,5)]=0.2", "maxsumW[c(1,4)]=0.9")
    #   portfolioConstraints(data, spec, constraints)

    # FUNCTION:

    # Already done ...
    if (class(constraints) == "fPFOLIOCON") return(constraints)

    # Missing target Return ...
    if (is.null(getTargetReturn(spec))) setTargetReturn(spec) <- NA

    # Handle NULL - A NULL  :
    if (is.null(constraints)) constraints="LongOnly"

    # Check Vector of Valid Strings - these are strings ...
    validStrings = c(
        "LongOnly", "Short",      # LongOnly and Short Notification
        "minW", "maxW",           # Box Constraints
        "minsumW", "maxsumW",     # left and right Sided Group Constraints
        "minB", "maxB",           # Covariance Risk Budgets
        "listF", "minF", "maxF",  # Nonlinear Functions Constraints
        "nCard", "minCard", "maxCard")
    if (any(constraints == "Short")) setSolver(spec) = "solveRshortExact"
    # usedStrings <- unique(sort(sub("\\[.*", "", constraints)))
    # checkStrings <- usedStrings %in% validStrings
    # check <- (sum(!checkStrings) == 0)
    # if (check) check <- "valid" else stop("Invalid Constraints String(s)")
    stringConstraints <- constraints
    # attr(stringConstraints, "control") = check

    # Data:
    Data <- portfolioData(data, spec)
    
    # Constraints:
    minW <- minWConstraints(Data, spec, constraints)
    maxW <- maxWConstraints(Data, spec, constraints)
    eqsumW <- eqsumWConstraints(Data, spec, constraints)
    minsumW <- minsumWConstraints(Data, spec, constraints)
    maxsumW <- maxsumWConstraints(Data, spec, constraints)
    minB <- minBConstraints(Data, spec, constraints)
    maxB <- maxBConstraints(Data, spec, constraints)
    listF <- listFConstraints(Data, spec, constraints)
    minF <- minFConstraints(Data, spec, constraints)
    maxF <- maxFConstraints(Data, spec, constraints)
    minBuyin <- minCardConstraints(Data, spec, constraints)
    maxBuyin <- maxCardConstraints(Data, spec, constraints)
    nCard <- nCardConstraints(Data, spec, constraints)
    minCard <- minCardConstraints(Data, spec, constraints)
    maxCard <- maxCardConstraints(Data, spec, constraints)

    if(is.null(minW)) minW = numeric()
    if(is.null(maxW)) maxW = numeric()
    if(is.null(eqsumW)) eqsumW = matrix(NA)
    if(is.null(minsumW)) minsumW = matrix(NA)
    if(is.null(maxsumW)) maxsumW = matrix(NA)
    if(is.null(minB)) minB = numeric()
    if(is.null(maxB)) maxB = numeric()
    if(is.null(maxsumW)) maxsumW = matrix(NA)
    if(is.null(minF)) minF = numeric()
    if(is.null(maxF)) maxF = numeric()
    if(is.null(minBuyin)) minBuyin = numeric()
    if(is.null(maxBuyin)) maxBuyin = numeric()
    if(is.null(nCard)) nCard = integer()
    if(is.null(minCard)) minCard = numeric()
    if(is.null(maxCard)) maxCard = numeric()
    
    # Return Value:
    new("fPFOLIOCON",
        stringConstraints=stringConstraints,
        minWConstraints=minW,
        maxWConstraints=maxW,
        eqsumWConstraints=eqsumW,
        minsumWConstraints=minsumW,
        maxsumWConstraints=maxsumW,
        minBConstraints=minB,
        maxBConstraints=maxB,
        listFConstraints=listF,
        minFConstraints=minF,
        maxFConstraints=maxF,
        minBuyinConstraints=minBuyin,
        maxBuyinConstraints=maxBuyin,
        nCardConstraints=nCard,
        minCardConstraints=minCard,
        maxCardConstraints=maxCard 
        )
}


################################################################################


minWConstraints <-
function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns a vector with min box constraints

    # Details:
    #   Takes care of "minW" strings, i.e. lower blounds
    #   W >= c

    # Arguments:
    #   data - a timeSeries or a fPFOLIODATA object
    #   spec - a fPFOLIOSPEC object
    #   constraints - a constraints string

    # Example:
    #   data <- as.timeSeries(data(LPP2005REC))[, 1:6]
    #   spec <- portfolioSpec()
    #   constraints <- c("minW[3:4]=0.1", "maxW[5:6]=0.8")
    #   minWConstraints(data, spec, constraints)

    # FUNCTION:

    # Settings:
    Data <- portfolioData(data, spec)
    if (class(data) == "fPFOLIODATA") data <- getSeries(Data)
    nAssets <- getNAssets(Data)
    assetsNames <- getUnits(Data)

    # Consider LongOnly:
    if("LongOnly" %in% constraints) {
        minW <- rep(0, nAssets)
        names(minW) <- assetsNames
        return(minW)
    }

    # Consider Unlimited Short:
    if("Short" %in% constraints) {
        minW <- rep(-Inf, nAssets)
        names(minW) <- assetsNames
        return(minW)
    }

    # Extract and Compose Vectors a_vec and b_vec:
    minW = rep(0, nAssets)
    names(minW) = assetsNames
    if (!is.null(constraints)) {
        nC <- length(constraints)
        what <- substr(constraints, 1, 4)
        for (i in 1:nC) {
            if (what[i] == "minW") eval(parse(text = constraints[i]))
        }
    }
    names(minW) <- assetsNames
    return(minW)

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


maxWConstraints <-
function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns a vector with max box constraints

    # Details:
    #   Takes care of "maxW" strings, i.e. upper bounds
    #   W >= c

    # Arguments:
    #   data - a timeSeries or a fPFOLIODATA object
    #   spec - a fPFOLIOSPEC object
    #   constraints - a constraints string

    # Example:
    #   data = as.timeSeries(data(LPP2005REC))[, 1:6]
    #   spec=portfolioSpec()
    #   constraints=c("minW[3:4]=0.1", "maxW[5:6]=0.8")
    #   maxWConstraints(data, spec, constraints)

    # FUNCTION:

    # Settings:
    Data <- portfolioData(data, spec)
    if (class(data) == "fPFOLIODATA") data <- getSeries(Data)
    nAssets <- getNAssets(Data)
    assetsNames <- getUnits(Data)

    # Consider LongOnly:
    if("LongOnly" %in% constraints) {
        maxW <- rep(1, nAssets)
        names(maxW) <- assetsNames
        return(maxW)
    }

    # Consider Unlimited Short:
    if("Short" %in% constraints) {
        maxW <- rep(Inf, nAssets)
        names(maxW) = assetsNames
        return(maxW)
    }

    # Extract and Compose Vectors a_vec and b_vec:
    maxW = rep(1, nAssets)
    names(maxW) <- assetsNames
    if (!is.null(constraints)) {
        nC <- length(constraints)
        what <- substr(constraints, 1, 4)
        for (i in 1:nC) {
            if (what[i] == "maxW") eval(parse(text = constraints[i]))
        }
    }
    names(maxW) = assetsNames
    return(maxW)

    # Return Value:
    invisible()
}


################################################################################


eqsumWConstraints <-
function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns a list with group equal matrix and vectors constraints

    # Details:
    #   Takes care of "eqsumW" strings
    #   A_eq W = c_eq

    # Arguments:
    #   data - a timeSeries or a fPFOLIODATA object
    #   spec - a fPFOLIOSPEC object
    #   constraints - a constraints string

    # Example:
    #   data = as.timeSeries(data(LPP2005REC))[, 1:6]
    #   spec=portfolioSpec(); setTargetReturn(spec) = mean(data)
    #   constraints="eqsumW[1:6]=1"
    #   eqsumWConstraints(data, spec, constraints)
    #   eqsumWConstraints(data, spec, constraints="LongOnly")
    #   eqsumWConstraints(data, spec, constraints=c("LongOnly","Partial"))

    # FUNCTION:

    # Get Statistics:
    Data <- portfolioData(data, spec)
    if (class(data) == "fPFOLIODATA") data <- getSeries(Data)
    targetReturn <- getTargetReturn(spec)[1]
    if (is.null(targetReturn)) {
        targetReturn = NA
        # stop("Target Return is Missing")
    }

    # Get Data:
    mu <- getMu(Data)
    nAssets <- getNAssets(Data)
    assetsNames <- getUnits(Data)

    # Target Return:
    Aeq <- matrix(mu, byrow = TRUE, ncol = nAssets)

    # Full or partial Investment?
    if ("partial" %in% tolower(constraints))
        fullInvest = FALSE else fullInvest = TRUE

    # Full Investment:
    #   - negative to handle better partial Investment in Rquadprog:
    if (fullInvest) Aeq <- rbind(Aeq, -rep(1, nAssets))

    # Dimension Names:
    colnames(Aeq) <- assetsNames
    if (fullInvest) {
        rownames(Aeq) <- c("Return", "Budget")
    } else {
        rownames(Aeq) <- "Return"
    }
    
    # RHS Vector:
    if (fullInvest) {
        ceq <- c(Return = targetReturn, Budget = -1)
    } else {
        ceq <- c(Return = targetReturn)
    }
    
    # Extract and Compose Matrix and Vector:
    what6 = substr(constraints, 1, 6)
    if (!is.null(constraints)) {
        nC = length(constraints)
        for (i in 1:nC) {
            if (what6[i] == "eqsumW")  {
                eqsumW = rep(0, times = nAssets)
                names(eqsumW) <- assetsNames
                eval(parse(text = constraints[i]))
                Aeq = rbind(Aeq, eqsumW = sign(eqsumW))
                a = strsplit(constraints[i], "=")[[1]][2]
                ceq = c(ceq, eqsumW = as.numeric(a))
            }
        }
    }
    eqsumW <- cbind(ceq, Aeq)

    # If target Return is missing:
    eqsumW <- na.omit(eqsumW)
    
    # Return Value:
    eqsumW
}


# ------------------------------------------------------------------------------


minsumWConstraints <-
function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns a list with group matrix and vectors constraints

    # Arguments:
    #   data - a timeSeries or a fPFOLIODATA object
    #   spec - a fPFOLIOSPEC object
    #   constraints - a constraints string$

    # Details:
    #   Takes care of "minsumW" strings
    #   a_vec <= A_mat W

    # Example:
    #   data = as.timeSeries(data(LPP2005REC))[, 1:6]
    #   spec=portfolioSpec(); setTargetReturn(spec) = mean(data)
    #   constraints=c("minsumW[2:3]=0.2", "minsumW[c(1,4:6)]=0.2")
    #   minsumWConstraints(data, spec, constraints)
    #   minsumWConstraints(data, spec)

    # FUNCTION:

    # Get Statistics:
    data <- portfolioData(data, spec)

    # Get Specifications:
    mu <- getMu(data)
    nAssets <- getNAssets(data)
    assetsNames <- getUnits(data)

    # Extrac and Compose Matrix and Vectors:
    what7 <- substr(constraints, 1, 7)

    if (!is.null(constraints)) {
        nC <- length(constraints)

        count <- 0
        Amat <- NULL
        avec <- NULL

        # Partial Investment:
        if ("partial" %in% tolower(constraints)) {
            Amat <- rbind(Amat, rep(1, times = nAssets))
            avec <- c(avec, 0)
        }

        for (i in 1:nC) {
            if (what7[i] == "minsumW")  {
                count = count + 1
                minsumW = rep(0, times = nAssets)
                names(minsumW) <- assetsNames
                eval(parse(text = constraints[i]))
                Amat = rbind(Amat, minsumW = sign(minsumW))
                a = strsplit(constraints[i], "=")[[1]][2]
                avec = c(avec, as.numeric(a))
            }
        }
        if (!is.null(Amat)){
            colnames(Amat) = assetsNames
            names(avec) = rep("lower", count)
        }

    }

    # Return Value:
    cbind(avec = avec, Amat = Amat)
}


# ------------------------------------------------------------------------------


maxsumWConstraints <-
function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns a list with group matrix and vectors constraints

    # Arguments:
    #   data - a timeSeries or a fPFOLIODATA object
    #   spec - a fPFOLIOSPEC object
    #   constraints - a constraints string$

    # Details:
    #   Takes care of "minsumW" and "maxsumW" strings
    #   a_vec <= A_mat W <= b_vec

    # Example:
    #   data = as.timeSeries(data(LPP2005REC))[, 1:6]
    #   spec=portfolioSpec(); setTargetReturn(spec) = mean(data)
    #   constraints=c("maxsumW[2:3]=0.7", "maxsumW[c(1,4:6)]=0.8")
    #   maxsumWConstraints(data, spec, constraints)
    #   maxsumWConstraints(data, spec)

    # FUNCTION:

    # Get Statistics:
    data <- portfolioData(data, spec)

    # Get Specifications:
    mu <- getMu(data)
    nAssets <- getNAssets(data)
    assetsNames <- getUnits(data)

    # Extract and Compose Matrix and Vectors:
    what7 = substr(constraints, 1, 7)

    if (!is.null(constraints)) {
        nC <- length(constraints)

        count <- 0
        Amat <- NULL
        avec <- NULL

        # Partial Investment:
        if ("partial" %in% tolower(constraints)) {
            Amat <- rbind(Amat, rep(1, times = nAssets))
            avec <- c(avec, 1)
        }

        for (i in 1:nC) {
            if (what7[i] == "maxsumW")  {
                count = count + 1
                maxsumW = rep(0, times = nAssets)
                names(maxsumW) <- assetsNames
                eval(parse(text = constraints[i]))
                Amat = rbind(Amat, maxsumW = sign(maxsumW))
                a = strsplit(constraints[i], "=")[[1]][2]
                avec = c(avec, as.numeric(a))
            }
        }
        if (!is.null(Amat)) {
            colnames(Amat) = assetsNames
            names(avec) = rep("upper", count)
        }
    }

    # Return Value:
    cbind(avec = avec, Amat = Amat)
}


################################################################################


minBConstraints <-
function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns a list with min risk budget constraints vectors

    # Arguments:
    #   constraints - a constraints string

    # Example:
    #   data = as.timeSeries(data(LPP2005REC))[, 1:6]
    #   spec=portfolioSpec()
    #   constraints=c("minB[3:4]=0.1","maxB[1:3]=0.3","maxB[c(4,6)]=0.4")
    #   minBConstraints(data, spec, constraints)
    #   minBConstraints(data, spec)

    # FUNCTION:

    # Create Data Object:
    Data <- portfolioData(data, spec)
    if (class(data) == "fPFOLIODATA") data <- getSeries(Data)

    # Get Specifications:
    nAssets <- getNAssets(Data)
    assetsNames <- getUnits(Data)

    # Extract and Compose Risk Budgets:
    minB = rep(-Inf, nAssets)
    names(minB) <- assetsNames
    if (!is.null(constraints)) {
        nC = length(constraints)
        what = substr(constraints, 1, 4)
        for (i in 1:nC) {
            if (what[i] == "minB") eval(parse(text = constraints[i]))
        }
    }

    # Return Value:
    minB
}


# ------------------------------------------------------------------------------


maxBConstraints <-
function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns a list with max risk budget constraints vectors

    # Arguments:
    #   constraints - a constraints string

    # Example:
    #   data = as.timeSeries(data(LPP2005REC))[,1:6]
    #   spec=portfolioSpec()
    #   constraints=c("minB[3:4]=0.1","maxB[1:3]=0.3","maxB[c(4,6)]=0.4")
    #   maxBConstraints(data, spec, constraints)
    #   maxBConstraints(data, spec)

    # FUNCTION:

    # Create Data Object:
    Data <- portfolioData(data, spec)
    if (class(data) == "fPFOLIODATA") data <- getSeries(Data)

    # Get Specifications:
    N <- nAssets <- getNAssets(Data)
    assetsNames <- getUnits(Data)

    # Extract and Compose Risk Budgets:
    maxB <- rep(1, N)
    names(maxB) <- assetsNames
    if (!is.null(constraints)) {
        nC = length(constraints)
        what = substr(constraints, 1, 4)
        for (i in 1:nC) {
            if (what[i] == "maxB") eval(parse(text = constraints[i]))
        }
    }

    # Return Value:
    maxB
}


################################################################################


listFConstraints <-
function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Nonlinear Constraints

    # Example: 
    #    maxdd <- function(x) max(drawdowns(x))
    #    listFConstraints(data <- NULL, constraints=c("minF=-0.04", "listF(maxdd)"))

    # FUNCTION:

    # Parse:
    nlin <- list()
    matched <- pmatch("listF" , constraints)
    if(!is.na(matched)) {
        constraints = paste("nlin = ", constraints[matched])
        constraints = sub("listF", "list", constraints)
        eval(parse(text = constraints))
    }
    
    # Return Value:
    return(nlin)
}


# ------------------------------------------------------------------------------


minFConstraints <-
function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Nonlinear Constraints

    # Example: 
    #    minFConstraints("minF=-0.04")

    # FUNCTION:

    # Parse:
    minF <- NULL
    matched <- pmatch("minF" , constraints)
    if(!is.na(matched)) eval(parse(text = constraints[matched]))
    
    # Return Value:
    return(minF)
}


# ------------------------------------------------------------------------------


maxFConstraints <-
function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Nonlinear Constraints

    # Example: 
    #    maxFConstraints(c("LongOnly", "maxF=0"))

    # FUNCTION:

    # Parse:
    maxF <- NULL
    matched <- pmatch("maxF" , constraints)
    if(!is.na(matched)) eval(parse(text = constraints[matched]))
    
    # Return Value:
    return(maxF)
}


################################################################################


minBuyinConstraints <-
function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Lower Buyin Constraints

    # Example: minBuyinConstraints(c("LongOnly", "minBuyin[1:3]=0.1")

    # FUNCTION:

    # Parse:
    Data <- portfolioData(data)
    nAssets <- getNAssets(Data)
    minBuyin <- rep(0, nAssets)
    matched <- pmatch("minBuyin", constraints)
    if(!is.na(matched)) eval(parse(text = constraints[matched]))
    names(minBuyin) <- getUnits(Data)
    
    # Return Value:
    return(minBuyin)
}


# ------------------------------------------------------------------------------


maxBuyinConstraints <-
function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Upper Buyininality Constraints

    # Example: maxBuyinConstraints(c("LongOnly", "maxBuyin[5]=0.9")

    # FUNCTION:

    # Parse:
    Data <- portfolioData(data)
    nAssets <- getNAssets(Data)
    maxBuyin <- rep(1, nAssets)
    matched <- pmatch("maxBuyin", constraints)
    if(!is.na(matched)) eval(parse(text = constraints[matched]))
    names(maxBuyin) <- getUnits(Data)
    
    # Return Value:
    return(maxBuyin)
}


################################################################################


nCardConstraints <-
function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Cardinality Constraints

    # Example: ncardConstraints(c("LongOnly", "ncard=4")

    # FUNCTION:

    # Parse:
    Data <- portfolioData(data)
    nAssets <- getNAssets(Data)
    nCard <- nAssets
    matched <- pmatch("nCard", constraints)
    if(!is.na(matched)) eval(parse(text = constraints[matched]))
    
    # Return Value:
    return(nCard)
}


# ------------------------------------------------------------------------------


minCardConstraints <-
function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Lower Cardinality Constraints

    # Example: minCardConstraints(c("LongOnly", "minCard[1:3]=0.1")

    # FUNCTION:

    # Parse:
    Data <- portfolioData(data)
    nAssets <- getNAssets(Data)
    minCard <- rep(0, nAssets)
    matched <- pmatch("minCard", constraints)
    if(!is.na(matched)) eval(parse(text = constraints[matched]))
    names(minCard) <- getUnits(Data)
    
    # Return Value:
    return(minCard)
}


# ------------------------------------------------------------------------------


maxCardConstraints <-
function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Upper Cardinality Constraints

    # Example: maxCardConstraints(c("LongOnly", "maxCard[5]=0.9")

    # FUNCTION:

    # Parse:
    Data <- portfolioData(data)
    nAssets <- getNAssets(Data)
    maxCard <- rep(1, nAssets)
    matched <- pmatch("maxCard", constraints)
    if(!is.na(matched)) eval(parse(text = constraints[matched]))
    names(maxCard) <- getUnits(Data)
    
    # Return Value:
    return(maxCard)
}


################################################################################
