"threshld" <-
function(x, t, hard = TRUE)
{
#
#  threshold the data x using threshold t
#  if hard=TRUE use hard thresholding
#  if hard=FALSE use soft thresholding
	if(hard) z <- x * (abs(x) >= t) else {
		z <- sign(x) * pmax(0, abs(x) - t)
	}
	return(z)
}
