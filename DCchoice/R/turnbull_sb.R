# Kaplan-Meier-Turnbull nonparametric approach to analyze 
# single-bounded dichotomous choice contingent valuation data

turnbull.sb <- function(formula, data, subset, conf.int = FALSE, B = 200, conf.level = 0.95, timeMessage = FALSE, ...){
    if(missing(data)) data <- environment(formula)
    
    if(length(formula[[2]]) != 1) stop("something is wrong with formula")
    
    mf <- match.call(expand.dots = TRUE)
    m <- match(c("formula", "data", "subset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$formula <- formula
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    original.data <- data
    data <- mf

    # removing observations with missing values
    na.num <- max(sum(as.numeric(is.na(data))))
    if(na.num != 0){ 
        d1 <- nrow(data)
        data <- na.omit(data)
        d2 <- nrow(data)
        warning(paste("Missing values detected.", d1 - d2, "rows are removed.", sep = " "))
    }
    
    # defining the dependent variable 
    lhs1 <- formula[[2]]      # extracting the name of the acceptance/rejection variable from the formula supplied
    y1 <- eval(lhs1, data)    # yes/no to the bids
    
    nobs <- length(y1)        # the number of observations
    
    P1 <- formula[[3]]        # retrieving the name of the bid variable from the formula supplied
    first <- eval(P1, data)   # the first stage bids
    
    # making dummy variables for the yes/no variable
    if(is.factor(y1)){   # when the yes/no variables are defined as factor
        y <- ifelse(y1 == "yes", 1, 0)
    } else {
        y <- y1
    }
    
    left <- ifelse(y == 1, first, 0)         # lower bound of WTP
    right <- ifelse(y == 1, Inf, first)      # upper bound of WTP
    unq.bid <- sort(unique(c(left, right)))  # unique bids including Inf
    
    turnbull <- icfit(L = left, R = right, conf.int  = conf.int, control = icfitControl(timeMessage = timeMessage, B = B, conf.level = conf.level))   # estimating nonparametric survival function. icfit function is defined in interval package
    
    # arranging outcomes into a single list variable
    output <- list(
        left = left,
        right = right,
        turnbull = turnbull,
        unq.bid = unq.bid
    )
    
    class(output) <- "turnbull"
    return(output)
    
}

