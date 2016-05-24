###############################################################################
#
# Function: plot_cs
#
# plot statistics
#
# Function calls: 	
#					
#
# Parameters:
#    	
#
# Return values:
#
# Author:			Niroshan Nadarajah <niroshan.nadarajah@uni-duesseldorf.de>
#
# Last modified:	10/11/14
#
###############################################################################

plot_avg_cs <- function(csstats, column) {
	
	plot(density(csstats@average[,column]), type='l', main="Average Tajima's D", sub="with observed values in red")
	abline(v=csstats@obsVal[,column], col=2, lty=3)
}