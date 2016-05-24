# $LastChangedDate: 2010-02-13 19:41:48 +0000 (Sat, 13 Feb 2010) $
# $Rev: 4745 $
#
# Unit test for "show" methods
#
# Author: Francisco
###############################################################################

test.show.COPPosterior <- function()
{
	x <- get(load( file.path(BLCOPOptions("unitTestPath"), "copexample.RData") ))
	
	checkEquals(capture.output(show(x$posterior)),  
			    c("Asset set:  SP,FTSE,CAC,DAX ", 
				  "Views used to generate this posterior: ", 
			      "[1] \"1*DAX~unif:(min=-0.02,max=0)\"", 
				  "Number of simulations: 1000 "
	             ))
}