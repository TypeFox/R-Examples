# Mallow's Cp

`Cp` <- 
function (object, ..., dispersion = NULL) {
    .cp <- function(x, disp) {
		rss <- deviance(x)
		df.r <- df.residual(x)
		if(is.null(df.r)) {
			ll <- logLik(x)
			df.r <- attr(ll, "nobs") - attr(ll, "df") + 1
		}
		scale <- if (!is.null(disp)) disp else if 
			(family(x)$family %in% c("poisson", "binomial")) 1 else if
			(df.r > 0) rss / df.r else NaN
		rss + 2 * scale  * (nobs(x) - df.r)
	}
    if (!missing(...)) {
        cps <- vapply(list(object, ...), .cp, 1L, disp = dispersion)
        val <- data.frame(Cp = cps)
        Call <- match.call()
        row.names(val) <- as.character(Call[-1L])
        val
    } else .cp(object, dispersion)
}