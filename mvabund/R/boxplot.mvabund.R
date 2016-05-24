################################################################################
# BOXPLOT.MVABUND: Display the abundances of several species in one plot or    #
# compare them for two objects                                                 #
################################################################################
boxplot.mvabund <- function (x, 
			y, 
			range = 1.5, 
			names=NULL,
    			at = NULL, 
			n.vars=min(12,NCOL(x)),
    			overall.main="Boxplot", 
			var.subset=NA, 
			transformation="log", 
			...)
{

	#Pass all arguments to the boxplot default call.
	default.boxplot.mvabund(x=x, 
			y=y, 
			range = range, 
#			names = names,
    			at = at, 
			n.vars = n.vars,
    			overall.main = overall.main, 
			var.subset = var.subset, 
			transformation = transformation, 
			...)

} 
	
