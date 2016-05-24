indiv.1m.ssc <- function(method, ES, cor, power = 0.8, alpha = 0.05, alternative = "two.sided", tol = 0.0001, maxiter = 1000, tol.uniroot = 0.0001) {

    if (missing(method)) stop("Missing 'method' argument.")
    if (missing(ES)) stop("Missing 'ES' argument.")
    if (missing(cor)) stop("Missing 'cor' argument.")
    if (class(method)!="character") stop("The 'method' argument should be of type character.")
    if (class(alternative) != "character") stop("The 'alternative' argument should be of type character.")	
    if(!is.matrix(cor)) stop("The 'cor' argument should be a matrix.")
    if(!is.vector(ES)) stop("The 'ES' argument should be a vector.")
  
  # Test of the function
    if( (alternative != "less") && (alternative != "greater") && (alternative != "two.sided") ) stop("The 'alternative' argument is misspecified.")
    if( (method != "Asympt") & (method != "Known") & (method != "Unknown") & (method != "UnKnown")) stop("The 'method' argument is misspecified.")	
    if ( (power < 0) || (power > 1)) stop("The 'power' argument should be between 0 and 1.")
    if ((alpha < 0) || (alpha > 1)) stop("The 'alpha' argument should be between 0 and 1.")
    if (!is.numeric(power)) stop("The 'power' argument should be numeric.")
    if (!is.numeric(alpha)) stop("The 'alpha' argument should be numeric.")  

    result <- .all.1m.ssc(method, ES, cor, power, alpha, alternative, tol, maxiter, tol.uniroot)
    
    return(result)

}

# See the paper 'Power and Sample Size Determination in Clinical Trials with
# Multiple Primary Continuous Correlated Endpoints', Journ. Biopharm. Stat., 2014, pages 384--385.
.all.1m.ssc <- function(method, ES, cor, power = 0.8, alpha = 0.05, alternative = "two.sided", tol = 0.0001, maxiter = 1000, tol.uniroot = 0.0001) {
    
    # 'ES':  Effect Size
    # 'cor': correlation matrix R
    # 'power': desired power 1 - \beta
    # 'alpha': FWER significance level
    
    # Initialization:
    alpha1 <- alpha
    stop.crit <- 1
    iter <- 0

    while ((stop.crit > tol) && (iter < maxiter)) {

        # Computation of the value of sample size n to obtain a power 'power', using individual critical values
        # given through 'lower' and 'upper' (namely, c{_alpha1 (or c_{\alpha1/2}):
        n <- .ncompute(method, ES, cor, power, alternative, alpha1)
        
        # Computation, for a given value of 'n', of the value alpha2 such that the value of FWER computed using 
        # the individual critical values qnorm(1 - alpha2) (which are noted c_{alpha1}) is equal to alpha:
        alpha2 <- .risk(method, n, alpha, length(ES), cor, alternative, tol.uniroot)

        # We will stop if alpha2 does not differ that much of alpha1:
        stop.crit <- abs(alpha1 - alpha2)
        iter <- iter + 1

        # Reinitialization of parameters
        alpha1 <- alpha2

    }

    return(list("Adjusted Type-I error rate" = alpha1, "Sample size" = n))
}


# Definition of power function
.powerfunc.indiv <- function(n, method = method, ES, cor, alternative, alpha1) {

    if ((method == "Known") || (method == "Asympt")) {
        myq <- qnorm
        mypm <- function(lower, upper) pmvnorm(lower, upper, mean = ES * sqrt(n / 2), sigma = cor)
    }
    if ((method == "Unknown") || (method == "UnKnown")) {
        df <- 2 * n - 2
        myq <- function(x) qt(x, df)
        mypm <- function(lower, upper) pmvt(lower, upper, df = df, delta = ES * sqrt(n / 2), sigma = cor)
    }
        
### !!!! Discuss with Jeremie and Ben!!! Jeremie had put two cases (for Unknown). If n<5 he used qnorm() instead of qt()
### What was the reason????        
        
    if (alternative == "two.sided") {
        lower <- -myq(1 - alpha1 / 2)
        upper <- myq(1 - alpha1 / 2)
    }
    if (alternative == "less") {
        lower <- -myq(1 - alpha1)
        upper <- Inf
    }
    if (alternative == "greater") {
        lower <- -Inf
        upper <- myq(1 - alpha1)
    }

    # Equation (7), page 382, of the power 1 - \beta in our 2014 paper with J. Riou and B. Liquet
    pow <- 1 - mypm(lower = rep(lower, length(ES)), upper = rep(upper, length(ES)))

    return(pow)
}


.f.n.indiv <- function(n, method = method, ES, cor, power, alternative, alpha1)
    return(.powerfunc.indiv(n, method = method, ES = ES, cor = cor, alternative, alpha1) - power)

# Definition of sample size determination function
.ncompute <- function(method, ES, cor, power, alternative, alpha1) {

    # Stop if sample size is larger than 10,000
    if (.f.n.indiv(10000, method = method, ES = ES, cor = cor, power = power, alternative = alternative, alpha1) < 0)
        stop("Sample size is larger than 10,000.")

    # Research of the sample size n
    n.value <- uniroot.integer(.f.n.indiv, interval = c(2, 10000), pos.side = TRUE, method = method, ES = ES, cor = cor,
                               power = power, alternative = alternative, alpha1 = alpha1)$root
    
    return(n.value)
}


# Definition of FWER function
.fwerfunc.indiv <- function(alpha, method, n, lenES, cor, alternative) {

    if ((method == "Known") || (method == "Asympt")) {
        myq <- qnorm
        mypm <- function(lower, upper) pmvnorm(lower, upper, mean = rep(0.0, lenES), sigma = cor)
    }
    if ((method == "Unknown") || (method == "UnKnown")) {
        df <- 2 * n - 2
        myq <- function(x) qt(x, df)
        mypm <- function(lower, upper) pmvt(lower, upper, df = df, delta = rep(0.0, lenES), sigma = cor)
    }

    if (alternative == "two.sided") {
        lower <- -myq(1 - alpha / 2)
        upper <- myq(1 - alpha / 2)
    }
    if (alternative == "less") {
        lower <- -myq(1 - alpha)
        upper <- Inf
    }
    if (alternative == "greater") {
        lower <- -Inf
        upper <- myq(1 - alpha)
    }

    # Equation (6), page 381, of the FWER in our 2014 paper with J. Riou and B. Liquet
    fwer <- 1 - mypm(lower = rep(lower, lenES), upper = rep(upper, lenES))
    
    return(fwer)
}

# Definition of the equation
.f.alpha <- function(alpha, method, alpha1, n, lenES, cor, alternative) return(.fwerfunc.indiv(alpha, method, n, lenES, cor, alternative) - alpha1)

# Function for Type-I error Determination
.risk <- function(method, n, alpha1, lenES, cor, alternative, tol.uniroot) {

    # Type I error rate computation
    alpha.value <- uniroot(.f.alpha, c(0, 1), method = method, alpha1 = alpha1, n = n, lenES = lenES,
                           cor = cor, alternative = alternative, tol = tol.uniroot)$root

    return(alpha.value)
    
}
