#################################################################
# 
# File:         ad.test.statistic.r
# Purpose:      Calculates the statistic for the Anderson Darling GoF test
#
# Created:      20090625
# Author:       Carlos J. Gil Bellosta
#
# Modifications: 
#
#################################################################

ad.test.statistic <- function( x )
{
    tmp <- x * ( 1 - rev( x ) )
    tmp <- ( 2 * seq(x) - 1 ) * log( tmp )
    tmp <- - mean( tmp ) - length( x )
}
