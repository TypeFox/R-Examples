bonferroni.1m.ssc <- function(mean.diff, sd, cor, power = 0.8, alpha = 0.05, alternative = "two.sided") {

    if (missing(mean.diff)) stop("Missing 'mean.diff' argument.")
    if (missing(sd)) stop("Missing 'sd' argument.")
    if (missing(cor)) stop("Missing 'cor' argument.")	
    if( (alternative != "less") && (alternative != "greater") && (alternative != "two.sided") )
        stop("The 'alternative' argument is misspecified.")
    if (class(alternative) != "character") stop("The 'alternative' argument must be of type character.")
    if (!is.vector(mean.diff)) stop("The 'mean.diff' argument should be a vector.")
    if (!is.vector(sd)) stop("The 'sd' argument should be a vector.")
    if (!is.matrix(cor)) stop("The 'cor' argument should be a matrix.")	
    if ((power < 0) || (power > 1)) stop("The 'power' argument should be between 0 and 1.")
    if ((alpha < 0) || (alpha > 1)) stop("The 'alpha' argument should be between 0 and 1.")
    if (!is.numeric(power)) stop("The 'power' argument should be numeric")
    if (!is.numeric(alpha)) stop("The 'alpha' argument should be numeric.")
		
    # Stop if sample size is upper than 10000 people
    if (.f.n.bonferroni(10000, power, alpha, sd, cor, mean.diff, alternative) < 0) stop("The sample size is larger than 10,000.")

    # Sample size computation
    n.value <- uniroot.integer(.f.n.bonferroni, c(2, 10000), pos.side = TRUE, power = power, alpha = alpha, sd = sd, cor = cor, mean.diff = mean.diff, alternative = alternative)$root
    return(n.value)
}

.power.bonferroni <- function(n, alpha, sd, cor, mean.diff, alternative) {
    
    # Covariance matrix
    matvar <- diag(sd) %*% cor %*% diag(sd)

    lenmd <- length(mean.diff)
    
    if (alternative == "two.sided") {
        lower <- -sd * qnorm(1 - alpha / (2 * lenmd), 0, 1) * sqrt(2 / n)
        upper <- sd * qnorm(1 - alpha / (2 * lenmd), 0, 1) * sqrt(2 / n)
    }
    if (alternative == "less") {
        lower <- -sd * qnorm(1 - alpha / lenmd, 0, 1) * sqrt(2 / n)
        upper <- rep(Inf, lenmd)
    }
    if (alternative == "greater") {
        lower <- rep(-Inf, lenmd)
        upper <- sd * qnorm(1 - alpha / lenmd, 0, 1) * sqrt(2 / n)
    }

                                        # Power Computation
    power <- 1 - pmvnorm(lower, upper, mean.diff, sigma = matvar * (2 / n))
    return(as.numeric(power))
}

.f.n.bonferroni <- function(n, power, alpha, sd, cor, mean.diff, alternative) {
    return(.power.bonferroni(n, alpha, sd, cor, mean.diff, alternative) - power)
}
