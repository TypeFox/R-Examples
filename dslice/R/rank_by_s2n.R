rank_by_s2n <- function(expmat, label)
{
	freq <- table(label)
	if((freq[1] < 2) || (freq[2] < 2)){
		stop("There is not enough sample for calculating singal to noise ratio")
	}
	x0 <- expmat[,which(label == 0)]
	x1 <- expmat[,which(label == 1)]
	m0 <- apply(x0, 1, mean)
	m1 <- apply(x1, 1, mean)
	sd0 <- apply(x0, 1, sd)
	sd1 <- apply(x1, 1, sd)
	s2n <- (m0 - m1) / (sd0 + sd1)
	grank <- rank(s2n)
	return(grank)
}