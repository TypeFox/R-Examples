summary.YPmodel <-
function(object, ...)
{
	YPmodelResult <- object

	Data <- YPmodelResult$Data
	Estimate <- YPmodelResult$Estimate
	IntervalBands <- YPmodelResult$IntervalBands
	LackFitTest <- YPmodelResult$LackFitTest
	Adlgrk <- YPmodelResult$Adlgrk

cat("Details of model of Yang and Prentice  \n")
summary.YPmodel.overall(object=Data)
summary.YPmodel.estimate(object=Estimate, interval=0)
summary.YPmodel.IntervalBands(IntervalBands)
summary.YPmodel.lackfittest(LackFitTest)
summary.YPmodel.adlgrk(Adlgrk)
cat("\n-------------------------------------------------------------------------------------------------------------  \n")

#plot.YPmodel(YPmodelResult)

}
