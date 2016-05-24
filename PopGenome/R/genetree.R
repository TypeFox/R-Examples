################################################
#
# Function: genetree
#
# This function uses Hudsons ms to generate samples under the neutral model
#
# Function calls: 	
#					progressBar
#
# Parameters:
#      model:   	
#
# Return values:
#
# Author:			Niroshan Nadarajah <niroshan.nadarajah@uni-duesseldorf.de>
#
# Last modified:	10/11/01
#
################################################

#genetree <- function(locus, n) {
	
#	if (class(locus)[1] != "locstats")
#		stop("First parameter must be an object of class 'locstats'")
#	
#	if ( !is.numeric(n) || n < 1)
#		stop("Please specify row number to print!")
	
#	library("ape")
#	phylo <- read.tree("", locus@trees[n])
#	plot(phylo, edge.color=c("black", "blue"))
#	title("Gene tree")
#}
