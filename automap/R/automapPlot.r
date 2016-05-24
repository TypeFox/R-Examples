automapPlot = function(plot_data, zcol, col.regions, ...)
# A fucntion to plot the results from autoKrige. Its purpose is to provide
# a custom colorscale usefull for the kriging results. 
{
	if(missing(zcol)) zcol = names(plot_data)
	if(missing(col.regions)) col.regions = c("#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#0C2C84","#5A005A")
	p = spplot(plot_data,
		zcol = zcol,
		col.regions=col.regions,
		cuts = length(col.regions) - 1,
		...)					

    return(p)
}