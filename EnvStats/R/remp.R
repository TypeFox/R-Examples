remp <-
function (n, obs) 
{
    ln <- length(n)
    if (ln < 1) 
        stop("'n' must be non-empty.")
    if (ln > 1) 
        n <- ln
    else {
        if (is.na(n) || n <= 0 || n != trunc(n)) 
            stop("'n' must be a positive integer or vector.")
    }
    if (!is.vector(obs, mode = "numeric")) 
        stop("'obs' must be a numeric vector")
    if ((bad.obs <- sum(!(obs.ok <- is.finite(obs)))) > 0) {
        is.not.finite.warning(obs)
        obs <- obs[obs.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'obs' removed."))
    }
    if (length(obs) < 2) 
        stop("'obs' must be a numeric vector with at least 2 non-missing elements")
    sample(obs, size = n, replace = TRUE)
}
