SPRT.default <-
function(distribution = "bernoulli", type1 = 0.05, type2 = 0.20, h0, h1, values = NULL, n = NULL, k = NULL) {

#SPRT <- function(distribution = "bernoulli", type1 = 0.05, type2 = 0.20, h0, h1, values = NULL, k = NULL, n = NULL) {
    
    ##################
    # Mandatory inputs
    
    # Distrtribution
    distribution <- as.character(distribution)
    
    if (length(distribution) != 1 || !(distribution %in% c("normal", "bernoulli", "poisson", "exponential"))) {
        stop("Set the distribution to one of: \"normal\", \"bernoulli\", \"poisson\", \"exponential\".")
    }
    
    # Errors
    type1 <- as.numeric(type1)
    type2 <- as.numeric(type2)
    
    (type1 > 0 & type1 < 1) || stop("Set the Type I error, \"type1\", between 0 and 1.")
    (type2 > 0 & type2 < 1) || stop("Set the Type II error, \"type2\", between 0 and 1.")
    
    # h0 and h1
    !is.null(h0) || stop("Please specify \"h0\", the parameter value associated with the null hypothesis.")
    !is.null(h1) || stop("Please specify \"h1\", the parameter value associated with the alternative hypothesis.")
    
    h0 <- as.numeric(h0)
    h1 <- as.numeric(h1)
    
    if (C.fn(distribution, h0) > C.fn(distribution, h1)) {
        stop("Please invert h0 and h1 to continue.")    
    }
    
    #################
    # Optional inputs
    
    # Random variable
    if (!is.null(values)) {
        
        # Remove missing values
        indices <- is.na(values)
        if (sum(indices) > 0) values <- values[!indices]
        rm(indices)
        
        # Number of observations
        n <- length(values)
        
        if (length(n) == 0) {
            values <- n <- k <- NULL
            warning("The random variable vector you supplied only contained NA values.")
            
        } else {
            # Sum the random variable
            k <- sum(values, na.rm = TRUE)
        }
        
        # Assess n and k, if supplied
    } else {
        
        is.null(n) || (n <- as.numeric(n))
        is.null(k) || (k <- as.numeric(k))
        
        # Avoid inversion
        if (is.numeric(n) & is.numeric(k)) {
            if (k > n & distribution == "bernoulli") {
                stop("\"k\" cannot be greater than \"n\" for a random variable with a Bernoulli distribution. \"k\" is the number of successes in \"n\" trials or observations.")
            }
        }
    }
    
    #######
    # LLR #
    #######
    
    # Wald boundaries
    boundary.A <- waldBoundary(type1=type1, type2=type2, boundary="A", log=TRUE)
    boundary.B <- waldBoundary(type1=type1, type2=type2, boundary="B", log=TRUE)
    
    # LLR coefficients
    llr.coefficients <- llr.fn(distribution=distribution, h0=h0, h1=h1)
    
    # LLR function
    closure.llr <- function(n.coefficient, k.coefficient) {
        function(n, k) {
            n * n.coefficient + k * k.coefficient
        }
    }
    
    llr.fn <- closure.llr(n.coefficient = llr.coefficients[["n"]], k.coefficient = llr.coefficients[["k"]])
    
    #############
    # SUCCESSES #
    #############
    
    # Success boundaries
    k.boundaries <- boundary.fn(distribution=distribution, type1=type1, type2=type2, h0=h0, h1=h1)
    
    # H0 and H1 functions
    closure.line <- function(intercept, slope) {
        function(n) {
            intercept + n * slope
        }
    }
    
    h0.fn <- closure.line(intercept = k.boundaries[1,1], slope = k.boundaries[1,2])
    h1.fn <- closure.line(intercept = k.boundaries[2,1], slope = k.boundaries[2,2])
    
    ########
    # SPRT #
    ########
    
    # Decision
    llr <- decision <- interpretation <- NULL
    
    if (is.numeric(n) & is.numeric(k)) {
        
        # Calculate the log likelihood ratio
        llr <- llr.fn(n=n, k=k)
        
        if (llr >= boundary.A) {
            decision <- FALSE
            interpretation <- "Accept H1"
            
        } else if (llr <= boundary.B) {
            decision <- TRUE
            interpretation <- "Accept H0"
            
        } else {
            decision <- NA
            interpretation <- "Continue testing"
        }
    }
    
    ##########
    # VALUES #
    ##########
    
    data.llr <- data.sum <- NULL
    
    if (!is.null(values)) {
        
        # LLR vs Wald boundaries
        data.llr <- data.frame(n = seq_along(values),
                               values = values,
                               k = cumsum(values),
                               wald.B = boundary.B,
                               wald.A = boundary.A)
        
        data.llr$llr <- llr.fn(n = data.llr$n, k = data.llr$k)
        
        # Successes vs H0 and H1 success boundaries
        data.sum <- data.frame(n = seq_along(values),
                               values = values,
                               k = cumsum(values))
        
        data.sum <- transform(data.sum,
                              h0 = h0.fn(n = n),
                              h1 = h1.fn(n = n))
    }
    
    ########
    # Output
    
    output <- list(distribution = distribution,
                   n = n,
                   k = k,
                   h0 = h0,
                   h1 = h1,                   
                   wald.A = boundary.A,
                   wald.B = boundary.B,
                   k.boundaries = k.boundaries,
                   llr.coefficients = llr.coefficients,
                   llr = llr,
                   decision = decision,
                   interpretation = interpretation,
                   data.llr = data.llr,
                   data.sum = data.sum,
                   llr.fn = llr.fn,
                   h0.fn = h0.fn,
                   h1.fn = h1.fn)
    
    class(output) <- "SPRT"
    
    output
}
