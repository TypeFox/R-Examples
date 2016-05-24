`find.ice` <-
function(aic)
{
	aic.min <- apply(aic, MARGIN = 1, FUN = min)
	aic.ind <- matrix(aic.min, nrow = nrow(aic), ncol = ncol(aic)) == aic
	aic.ind <- aic.ind & 1 == t(apply(aic.ind, MARGIN = 1, cumsum))
	ice <- t(col(aic))[t(aic.ind)] - 1
	ice
}

