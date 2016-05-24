eqsCorr <-
function(eqs)
{
	options(warn=-1)
	eqs		<- eqsCov(eqs)
	eqs		<- cov2cor(eqs)
	return(eqs)
}
