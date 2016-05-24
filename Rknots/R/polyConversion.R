##Script generated in:
# 2011
# 2:22:13 PM
#by: 
# Author: Federico Comoglio @ D-BSSE, ETH Zurich
###############################################################################

HOMFLY2Alexander <- function(HOMFLY) {
	Var('l')
	Var('m')
	Var('x')
	HOMFLY <- parseToSympy(HOMFLY)
	toeval <- paste("((", HOMFLY, ")", 
			".subs(l, 1)).subs(m, (sqrt(x) - 1 / sqrt(x))).expand()", sep = "")
	return( parseToR( sympy(toeval) ))
}

HOMFLY2Jones <- function(HOMFLY) {
	Var('l')
	Var('m')
	Var('t')
	HOMFLY <- parseToSympy(HOMFLY)
	toeval <- paste("((", HOMFLY, ")", 
			".subs(l, t**-1)).subs(m, (sqrt(t) - 1 / sqrt(t))).expand()", sep = "")
	return( parseToR( sympy(toeval) ))
}

HOMFLY2mirror <- function(HOMFLY) {
	Var('l')
	Var('m')
	HOMFLY <- parseToSympy(HOMFLY)
	toeval <- paste("((", HOMFLY, ")", 
			".subs(l, -l**-1)).expand()", sep = "")
	return( parseToR( sympy(toeval) ))
}


