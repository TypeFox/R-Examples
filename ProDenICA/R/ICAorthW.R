ICAorthW <-
function(W)
{
	sW <- svd(W)
	sW$u %*% t(sW$v)
}

