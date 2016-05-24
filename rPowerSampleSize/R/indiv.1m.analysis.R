indiv.1m.analysis <- function(method, XC, XT, varX = NULL, alpha = 0.05, alternative = "two.sided", n = NULL) {

# EST - CE QUE n n'est pas egal a nrow(XC) ???
    
    if (missing(method)) stop("Missing 'method' argument.")
    if (missing(XC)) stop("Missing 'XC' argument.")
    if (missing(XT)) stop("Missing 'XT' argument.")
    if (class(method) != "character") stop("The 'method' argument should be of type character.")
    if( (method != "Asympt") && (method != "Known") && (method != "UnKnown"))
        stop("The 'method' argument is misspecified.")
    if( (alternative != "less") && (alternative != "greater") && (alternative != "two.sided") )
        stop("The 'alternative' argument is misspecified.")
    if (class(alternative) != "character") stop("The 'alternative' argument should be of type a character.")
    if ((alpha < 0) || (alpha > 1)) stop("The 'alpha' argument should be between 0 and 1.")
    if (!is.numeric(alpha)) stop("The 'alpha' argument should be numeric.")
    n <- nrow(XC)
    if (nrow(XT) != n) stop("'XC' and 'XT' should have the same number of rows.")
    
    # Computation of the correlation matrix
    cor <- cor(rbind(XC, XT))
	
    # Computation of the test statistics, and of the corresponding unadjusted p-values
    Indivs <- .test.indiv(method, XC, XT, n, varX, alternative)
    
    nb.stat <- length(Indivs$pval)
    m <- ncol(XT)
    
    p.val <- rep(NA, nb.stat)
    for (i in 1:nb.stat) {
        p.val[i] <- .fwer.pval.indiv(as.numeric(Indivs$stat[i]), m, alternative, method, cor, n)
    }
    
    names(p.val) <- 1:ncol(XC)
    
    result <- list("UnAdjPvalue" = Indivs$pval, "AdjPvalue" = p.val)
    class(result) <- c("indiv.1m.analysis")
    return(result)
    
}

# Function for the computation of the adjusted p-values
.fwer.pval.indiv <- function(stat.obs, m, alternative, method, cor, n) { 

    if (alternative == "two.sided") {
        a <- 1 ; b <- -1
        lower <- rep(-abs(stat.obs), m)
        upper <- rep(abs(stat.obs), m)
    }
    if (alternative == "less") {
        a <- 0 ; b <- 1
        lower <- rep(-Inf, m)
        upper <- rep(stat.obs, m)
    }
    if (alternative == "greater") {
        a <- 0 ; b <- 1
        lower <- rep(stat.obs, m)
        upper <- rep(Inf, m)
    }

    if ((method == "Asympt") || (method == "Known")) {  
        fwer <- a + b * pmvnorm(lower = lower, upper = upper, mean = rep(0.0, m), sigma = cor)
    }
    if ((method == "UnKnown") || (method == "Unknown") ) { 
        fwer <- a + b * pmvt(lower = lower, upper = upper, delta = rep(0.0, m), sigma = cor, df = 2 * n - 2)
    }
    
    return(as.numeric(fwer))
}

.test.indiv <- function(method, XC, XT, n, varX = NULL, alternative = "two.sided") {

    if ((method == "Known") && is.null(varX)) stop("When 'method' = 'Known', you should provide the value of 'varX'")
    
    # Computation of the mean difference for each endpoint
    num <- colMeans(XT) - colMeans(XC)

    # Computation of the test statistic, and of the corresponding p-values
    if ((method == "Unknown") || (method == "UnKnown")) {
        Y <- rbind(XC, XT)
        varY <- apply(Y, MARGIN = 2, FUN = var)
        # Equation (5) page 381 in our 2014 paper:
        Stat <- num / sqrt((2 / n) * varY)

#### For the paper published in 2014, there was a mistake in the code below. The degree of freedom was incorrectly set to df = n.
        if (alternative == "two.sided") pval <- 2 * (1 - pt(abs(Stat), df = 2 * n - 2))
        if (alternative == "greater") pval <- 1 - pt(Stat, df = 2 * n - 2)
        if (alternative == "less") pval <- pt(Stat, df = 2 * n - 2)
    }
    if (method == "Known") {
        # Equation (4) page 381 in our 2014 paper:
        Stat <- num / sqrt((2 / n) * varX)
        if (alternative == "two.sided") pval <- 2 * (1 - pnorm(abs(Stat)))
        if (alternative == "greater") pval <- 1 - pnorm(Stat)
        if (alternative == "less") pval <- pnorm(Stat)
    }
    if (method == "Asympt") {
        # Computation of the variances:
        varXC <- apply(XC, MARGIN = 2, FUN = var)
        varXT <- apply(XT, MARGIN = 2, FUN = var)
        # Equation (8) page 383 in our 2014 paper:
        Stat <- num / sqrt((varXT + varXC) / n)
        if (alternative == "two.sided") pval <- 2 * (1 - pnorm(abs(Stat)))
        if (alternative == "greater") pval <- 1 - pnorm(Stat)
        if (alternative == "less") pval <- pnorm(Stat)
    }
    
    result <- list(pval = pval, stat = Stat)
    return(result)
}
