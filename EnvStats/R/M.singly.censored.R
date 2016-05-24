M.singly.censored <-
function (data, left.censor, t.df = 3) 
{
    if ((bad.obs <- sum(!(ok <- is.finite(data)))) > 0) {
        is.not.finite.warning(data)
        data <- data[ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'data' removed."))
    }
    if (!is.vector(data, mode = "numeric")) 
        stop("'data' must be a numeric vector")
    if (length(data) < 2) 
        stop("'data' must have 2 or more observations")
    if (!(length(left.censor) == 1)) 
        stop("'left.censor' must be a constant")
    if (t.df < 1) 
        stop("'t.df' should be greater than or equal to 1")
    data.name <- deparse(substitute(data))
    censoring.name <- deparse(substitute(censored))
    censoring.side <- "Left"
    output <- fix0205(data, left.censor, alpha = t.df/2, beta = 1)
    ret.list <- list(distribution = "Normal", sample.size = output$sample[1], 
        censoring.side = censoring.side, censoring.levels = output$sample[5], 
        percent.censored = (100 * output$sample[3])/output$sample[1], 
        parameters = output$consistent.estimates, n.param.est = 2, 
        method = paste("M-estimator: t.df =", t.df), data.name = data.name, 
        censoring.name = censoring.name, bad.obs = bad.obs, var.cov.params = output$covariance.estimates)
    oldClass(ret.list) <- "estimateCensored"
    ret.list
}
