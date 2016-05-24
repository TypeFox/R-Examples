funopare.kernel <-
function(Response, CURVES, PRED, bandwidth, ..., kind.of.kernel = "quadratic", semimetric = "deriv")
{
    Response <- as.vector(Response)
	if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
	testfordim <- sum(dim(CURVES)==dim(PRED))==2
	twodatasets <- T
	if(testfordim) twodatasets <- sum(CURVES==PRED)!=prod(dim(CURVES))
	sm <- get(paste("semimetric.", semimetric, sep = ""))
	if(semimetric == "mplsr")
		SEMIMETRIC1 <- sm(Response, CURVES, CURVES, ...)
	else SEMIMETRIC1 <- sm(CURVES, CURVES, ...)
	kernel <- get(kind.of.kernel)
	KERNEL1 <- kernel(SEMIMETRIC1/bandwidth)
	KERNEL1[KERNEL1 < 0] <- 0
	KERNEL1[KERNEL1 > 1] <- 0
       diag(KERNEL1) <- 0
	RESPKERNEL1 <- KERNEL1 * Response
	Denom1 <- apply(KERNEL1, 2, sum)
	if(sum(Denom1 == 0) > 0)
	{
		return(list(Mse = 0))
	}
	else
	{
		NWweit = KERNEL1/Denom1;
		Response.estimated <- apply(RESPKERNEL1, 2, sum)/Denom1
		Mse.estimated <- sum((Response.estimated - Response)^2)/length(Response)
	}
	if(twodatasets) 
	{
		if(semimetric == "mplsr")
			SEMIMETRIC2 <- sm(Response, CURVES, PRED, ...)
		else SEMIMETRIC2 <- sm(CURVES, PRED, ...)
		KERNEL2 <- kernel(SEMIMETRIC2/bandwidth)
		KERNEL2[KERNEL2 < 0] <- 0
		KERNEL2[KERNEL2 > 1] <- 0
		Denom2 <- apply(KERNEL2, 2, sum)
		if(sum(Denom2 == 0) > 0)
		{
			return(0)
		}
		else
		{
			RESPKERNEL2 <- KERNEL2 * Response
			Response.predicted <- apply(RESPKERNEL2, 2, sum)/Denom2
			return(list(NWweit = NWweit, Estimated.values = Response.estimated, 
				Predicted.values = Response.predicted, band = bandwidth,
				Mse = Mse.estimated))
		}
	}
	else 
	{
		return(list(NWweit = NWweit, Estimated.values = Response.estimated, band = 
					bandwidth, Mse = Mse.estimated))
	}
}
