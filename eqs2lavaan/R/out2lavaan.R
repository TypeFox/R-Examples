out2lavaan		<- function(eqs)
{
	options(warn=-1)
	eqs		<- readLines(eqs, n=-1L)
	covi	<- eqsCov(eqs)
	desc	<- eqsDesc(eqs)
	return(list(covi,desc))
}