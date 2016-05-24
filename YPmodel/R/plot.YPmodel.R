plot.YPmodel <-
function(x, ...)
{

	YPmodelResult <- x

	Data <- YPmodelResult$Data
	IntervalBands <- YPmodelResult$IntervalBands
	LackFitTest <- YPmodelResult$LackFitTest
	Estimate <- YPmodelResult$Estimate

	par(mfrow=c(2,2))
	plot.YPmodel.survivor(Estimate)
	plot.YPmodel.IntervalBands(IntervalBands)
	plot.YPmodel.survf(LackFitTest)
	plot.YPmodel.martint(LackFitTest)

}
