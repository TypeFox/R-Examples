#---------------------------------------------------------------------------
#
#   This file holds the S4 class union definitions that are used in the 
#   sampSurf package classes.
#
#   classUnion...
#   1. numericNull
#
#Author...									Date: 26-Aug-2010
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
# used for solidType in downLog object...
#
setClassUnion('numericNULL', c('numeric', 'NULL'))

#
# general data frame class that allows 'missing' in function arguments... 
#
#setClassUnion('missDF', c('data.frame','missing'))

