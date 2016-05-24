################################################################################
# PLOTMVAFACTOR                                                                #
################################################################################
plotMvaFactor <- function(x, 
			y, 
			type="p", 
			main="Abundance", 
			n.vars= min(12,NCOL(x)),
			transformation="log", 
			legend = TRUE, ... ) 
{
	#Send default MVAFACTOR plot.
	default.plotMvaFactor(x, y, 
			type=type, 
			main=main, 
			n.vars=n.vars,
			transformation=transformation, 
			legend=legend, ...)
 
}

