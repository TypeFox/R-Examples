# Copyright (C) 2011 Jelmer Ypma. All Rights Reserved.
# This code is published under the GPL.
#
# File:   readASCGrid.R
# Author: Jelmer Ypma
# Date:   7 November 2011
#
# Description: function to read an .asc file with a sparse grid
#			   works for .asc files from www.sparse-grids.de
#
# Input: 
#     filename		name of the .asc file that you want to read
#					the extension .asc should be included
#	  dimension 	dimension of the grid that is being read
#
# Output: 
#     list with two elements:
#		nodes    	= matrix of nodes with dim columns 
#     	weights    = row vector of corresponding weights
#

readASCGrid <- function( filename, dimension ) {
	int.grid		<- matrix( scan( filename, sep=',' ), ncol=dimension+1, byrow=TRUE )
	nodes			<- matrix( int.grid[,-ncol(int.grid) ], ncol=ncol(int.grid)-1 )		# all columns expect last one
	weights			<- int.grid[, ncol(int.grid) ]										# last column
	
	return( list( "nodes" = nodes, "weights" = weights ) )
}
