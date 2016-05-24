################################################################################
# PLOT.MVABUND: Plot functions for mvabund objects (multivariate abundance data)
################################################################################
plot.mvabund <- function(x, 
				y, 
				type="p", 
				overall.main="",
				n.vars= 12,
				var.subset=NA,
  				transformation="log", 
				...) 
{ 
	#Pipe function call to the default plot call with many different options.
	default.plot.mvabund(x=x, 
				y=y, 
				type=type, 
				overall.main=overall.main,
				var.subset=var.subset,
				n.vars=n.vars,
  				transformation=transformation, 
				...)
}
###### END FUNCTION CALL ######
setMethod("plot", signature("mvabund", "mvabund"), plot.mvabund)
setMethod("plot", signature("mvabund", "missing"), plot.mvabund)
setMethod("plot", signature(x="ANY", y="mvabund"), plot.mvabund)
setMethod("plot", signature(x="mvabund",y="ANY"), plot.mvabund)

