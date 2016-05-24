plot.YPmodel.lackfittest <-
function(x=c(), Internal=c(), ...)
#######################################################
# Sep 20, 2015 - v0.1 Create
#######################################################
{
	#dev.new()
	#attach(mtcars)
	par(mfrow=c(1,2))
	plot.YPmodel.martint(x, Internal)
	plot.YPmodel.survf(x, Internal)
	#dev.off()
}
