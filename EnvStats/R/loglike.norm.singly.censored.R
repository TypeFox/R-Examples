loglike.norm.singly.censored <-
function (x, censored, censoring.side, mean, sd) 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector")
    if (!(is.vector(censored, mode = "numeric") || is.vector(censored, 
        mode = "logical"))) 
        stop("'censored' must be a logical or numeric vector")
    if (length(censored) != length(x)) 
        stop("'censored' must be the same length as 'x'")
    if (any(is.na(censored))) 
        stop("'censored' cannot contain missing values")
    if (is.numeric(censored)) {
        if (!all(censored == 0 | censored == 1)) 
            stop(paste("When 'censored' is a numeric vector, all values of", 
                "'censored' must be 0 (not censored) or 1 (censored)."))
        censored <- as.logical(censored)
    }
    n.cen <- sum(censored)
    if (n.cen == 0) 
        stop("No censored values indicated by 'censored'.")
    censoring.side <- match.arg(censoring.side, c("left", "right"))
    if (!is.vector(mean, mode = "numeric") || length(mean) != 
        1 || !is.finite(mean)) 
        stop("'mean' must be a non-missing, finite numeric scalar")
    if (!is.vector(sd, mode = "numeric") || length(sd) != 1 || 
        !is.finite(sd) || sd < 0) 
        stop("'sd' must be a non-missing, finite, positive numeric scalar")
    data.name <- deparse(substitute(x))
    censoring.name <- deparse(substitute(censored))
    if (any(is.na(x))) 
        statistic <- NA
    else {
        x.no.cen <- x[!censored]
        if (length(unique(x.no.cen)) < 2) 
            stop(paste("'x' must contain at least 2 non-missing,", 
                "uncensored, distinct values."))
        T1 <- unique(x[censored])
        if (length(T1) > 1) 
            stop(paste("More than one censoring level. ", "Use 'loglike.norm.multiply.censored'."))
        if (censoring.side == "left") {
            if (T1 > min(x.no.cen)) 
                stop(paste("For singly left-censored data,", 
                  "all uncensored observations must be bigger than", 
                  "or equal to the censoring level. ", "Use 'loglike.norm.multiply.censored'."))
        }
        else {
            if (T1 < max(x.no.cen)) 
                stop(paste("For singly right-censored data,", 
                  "all uncensored observations must be less than", 
                  "or equal to the censoring level. ", "Use 'loglike.norm.multiply.censored'."))
        }
        N <- length(x)
        n <- length(x.no.cen)
        stat1 <- lchoose(N, n.cen)
        con <- ifelse(censoring.side == "left", T1, -T1)
        stat2 <- n.cen * log(pnorm(con, mean = mean, sd = sd))
        stat3 <- (-n/2) * log(2 * pi) - n * log(sd) - (1/(2 * 
            sd^2)) * sum((x.no.cen - mean)^2)
        statistic <- stat1 + stat2 + stat3
    }
    names(statistic) <- "Log-Likelihood"
    parameters <- c(mean, sd)
    names(parameters) <- c("mean", "sd")
    list(distribution = "Normal", dist.abb = "norm", distribution.parameters = parameters, 
        statistic = statistic, sample.size = N, censoring.side = censoring.side, 
        censoring.levels = T1, percent.censored = (100 * n.cen)/N, 
        data.name = data.name, censoring.name = censoring.name)
}
